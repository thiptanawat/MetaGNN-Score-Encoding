#!/usr/bin/env python3
"""Build the independent-evaluation inputs from the public Copeland et al. 2023 sources.

Reads only public files: the authors' research compendium and the companion RNA-seq data
package, both pinned by commit. Produces the reaction-evidence matrix, the measured exchange
rates and the measured growth rates, each with source hashes.

Preprocessing interface, fixed here and recorded in the plan before any agreement is computed:

  * Counts are converted to counts per million within each library. Gene lengths are not
    published with these data, so no length normalization is applied and the result is NOT
    called TPM. The original study's GCRMA scale constant is not transferred; a new scale
    constant is fitted from these data alone, by the same rule the original used.
  * Recon3D gene identifiers carry a transcript suffix; the Entrez base before the first dot
    is the join key, matching the study's own model gene map.
  * Several rows may share an Entrez base. They are aggregated by maximum, which is the rule
    the original analysis used for probes sharing an Entrez identifier.
  * Gene-protein-reaction rules are evaluated with the study's own conservative parser, so an
    AND whose operand is unmeasured leaves the reaction unsupported rather than guessing.

No measured glucose or lactate value is read by this script.
"""
import argparse, csv, hashlib, json, os, re, subprocess, sys
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit(repo):
    return subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"]).decode().strip()


def read_rda(path):
    import rdata
    return rdata.conversion.convert(rdata.parser.parse_file(path))


def cat_values(x):
    """R factors arrive as pandas Categorical; plain vectors arrive as arrays."""
    return [str(v) for v in list(x)]


def load_counts(rnaseq_repo):
    """Raw integer counts, genes by samples, with the sample design table."""
    se = read_rda(Path(rnaseq_repo) / "data" / "lf_hyp_bay_rnaseq.rda")["lf_hyp_bay_rnaseq"]
    counts = np.asarray(se.assays.data.listData["counts"])
    genes = [str(g) for g in np.asarray(se.NAMES)]
    col = se.colData.listData
    samples = [str(s) for s in np.asarray(col["id"])]
    experiment = [str(v) for v in np.asarray(col["experiment"])]
    oxygen = cat_values(col["oxygen"])
    treatment = cat_values(col["treatment"])
    meta = {s: {"experiment": experiment[i], "oxygen": oxygen[i], "treatment": treatment[i]}
            for i, s in enumerate(samples)}
    return counts, genes, samples, meta


def load_growth_rates(compendium):
    """Measured per-condition growth rates for the lf_05-bay experiment, in per hour.

    R/data.R documents mu as "growth rate per hour" and X0 as "cell count at time 0". Growth
    rate is a declared model input here, not one of the predicted targets.
    """
    df = read_rda(Path(compendium) / "data" / "growth_rates.rda")["growth_rates"]
    df.columns = [str(c) for c in df.columns]
    keep = df[(df["cell_type"].astype(str) == "lf") &
              (df["experiment"].astype(str) == "05-bay")].copy()
    if len(keep) != 16:
        raise ValueError("Expected 16 lf/05-bay growth-rate records, found %d" % len(keep))
    keep = keep.astype({"X0": float, "mu": float})
    for c in keep.columns:
        if c not in ("X0", "mu"):
            keep[c] = keep[c].astype(str)
    return keep.to_dict("records")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study-root", type=Path, required=True)     # the pinned cpu_study directory
    ap.add_argument("--compendium", type=Path, required=True)
    ap.add_argument("--rnaseq", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(a.study_root))
    import cobra
    from src.gpr import evaluate

    counts, genes, samples, meta = load_counts(a.rnaseq)
    growth = load_growth_rates(a.compendium)
    if counts.shape != (len(genes), len(samples)):
        counts = counts.T
    if counts.shape[1] != 16:
        raise ValueError("Expected 16 RNA samples, found %d" % counts.shape[1])
    if not np.all(np.isfinite(counts)) or (counts < 0).any():
        raise ValueError("Counts must be finite and nonnegative")

    libsize = counts.sum(axis=0)
    cpm = counts / libsize * 1e6                      # outcome-blind, no gene-length assumption

    # aggregate to Entrez base by maximum, the rule the original analysis used for shared probes
    by_gene = {}
    for i, g in enumerate(genes):
        base = str(g).split(".")[0]
        if not base.isdigit():
            continue
        by_gene[base] = np.fmax(by_gene[base], cpm[i]) if base in by_gene else cpm[i].copy()

    model = cobra.io.load_json_model(str(a.study_root / "data/raw/Recon3D.json"))
    # reuse the study's own model-gene to Entrez map rather than re-deriving the identifier form
    gene_map = {r["model_gene_id"]: r["entrez_id"] for r in
                csv.DictReader(open(a.study_root / "manifests/model_gene_map.tsv"), delimiter="\t")}
    model_gene_ids = {g.id for g in model.genes}
    if not model_gene_ids <= set(gene_map):
        raise ValueError("Model gene identifiers are not covered by the study gene map")
    model_genes = {gid: gene_map[gid] for gid in model_gene_ids}
    measured = {gid: by_gene[base] for gid, base in model_genes.items()
                if base and base in by_gene}
    if not measured:
        raise ValueError("No model gene matched the expression matrix; check the identifier form")

    reaction_ids, rows, unsupported = [], [], []
    for r in model.reactions:
        if r.boundary:
            continue
        val = evaluate(r.gene_reaction_rule, measured, policy="conservative")
        if val is None:
            unsupported.append(r.id)
            continue
        reaction_ids.append(r.id); rows.append(np.asarray(val, dtype=float))
    A = np.vstack(rows)
    if (A < 0).any() or not np.isfinite(A).all():
        raise ValueError("Reaction evidence must be finite and nonnegative")

    k = float(np.median(A[A > 0]))                    # same rule as the original scale constant
    np.savez_compressed(a.out / "scores_copeland.npz", A=A, rxn=np.array(reaction_ids),
                        contexts=np.array(samples), cols=np.array(samples))

    internal = sum(1 for r in model.reactions if not r.boundary)
    with_rule = sum(1 for r in model.reactions if not r.boundary and r.gene_reaction_rule.strip())
    summary = {
        "source": {
            "compendium_commit": git_commit(a.compendium),
            "rnaseq_commit": git_commit(a.rnaseq),
            "compendium_url": "https://github.com/oldhamlab/Copeland.2023.hypoxia.flux",
            "rnaseq_url": "https://github.com/wmoldham/rnaseq.lf.hypoxia.molidustat",
            "article_doi": "10.7554/eLife.82597",
            "sra_bioproject": "PRJNA721596",
            "rnaseq_rda_sha256": sha256(Path(a.rnaseq) / "data" / "lf_hyp_bay_rnaseq.rda"),
            "recon3d_sha256": sha256(a.study_root / "data/raw/Recon3D.json"),
        },
        "expression": {
            "genes_in_matrix": len(genes), "samples": len(samples),
            "library_sizes_min": float(libsize.min()), "library_sizes_max": float(libsize.max()),
            "normalization": "counts per million within library; no gene-length normalization; not TPM",
            "entrez_bases_measured": len(by_gene),
            "model_genes_total": len(model_genes), "model_genes_measured": len(measured),
        },
        "reactions": {
            "internal_reactions": internal, "internal_with_gpr": with_rule,
            "supported_conservative": len(reaction_ids), "unsupported": len(unsupported),
            "policy": "conservative: an AND with an unmeasured operand leaves the reaction unsupported",
        },
        "scale_constant_k": k,
        "scale_constant_rule": "median of positive reaction evidence over all 16 samples of this dataset; the original GCRMA-derived constant is not transferred",
        "samples": meta,
        "growth_rate_records": len(growth),
        "growth_rate_units": "per hour, as documented in the compendium R/data.R",
        "no_target_outcome_read": True,
        "note": "Measured growth rate is a declared model input. Measured glucose and lactate rates are the prediction targets and are not read by this script.",
    }
    json.dump(summary, open(a.out / "prepare_summary.json", "w"), indent=1)
    with open(a.out / "growth_rates.tsv", "w", newline="") as fh:
        fields = sorted({k for r in growth for k in r})
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(growth)
    with open(a.out / "sample_table.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "experiment", "oxygen", "treatment", "library_size"])
        for i, s in enumerate(samples):
            w.writerow([s, meta[s]["experiment"], meta[s]["oxygen"], meta[s]["treatment"],
                        int(libsize[i])])
    print(json.dumps({kk: vv for kk, vv in summary.items() if kk != "samples"}, indent=1))


if __name__ == "__main__":
    main()
