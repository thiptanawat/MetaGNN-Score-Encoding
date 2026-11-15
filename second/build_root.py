#!/usr/bin/env python3
"""Second reconstruction: build a study root in which the pipeline of the main analysis runs on Human-GEM.

    python3 second/build_root.py --main MAIN_ROOT --humangem HUMAN_GEM_CHECKOUT --out SECOND_ROOT

The pipeline of the main analysis reads one model file, data/raw/Recon3D.json, opens its medium by the
identifiers EX_<stem>_e, and reads reaction evidence from data/scores_conservative.npz. This script
writes a second root with the same layout whose model is Human-GEM, so that the scripts under src/,
mechanism/ and closure/ run in it unchanged. Nothing here reads a measured exchange (no file under
data/ of the main root other than gene_expression.npz is opened).

What is changed, and only this:

1. Model. Human-GEM is read from its SBML file and written as COBRA JSON under the file name the
   pipeline reads. Stoichiometry, bounds, gene rules and the biomass reaction are untouched. Boundary
   reactions are renamed by Human-GEM's own cross-reference to BiGG: EX_<metBiGGID>_e, where metBiGGID
   is the BiGG identifier that model/metabolites.tsv gives for the exchanged metabolite, which is how
   the main model names its exchanges; where the metabolite has none, the reaction's own BiGG
   identifier (model/reactions.tsv, column rxnBiGGID). A boundary reaction takes that identifier when
   it is the only boundary reaction that would, and EX_<its Human-GEM id> otherwise. second/target_map.tsv
   also records, for every target, whether the reaction-level cross-reference gives the same identifier. The cobamide-coenzyme exchange keeps its Human-GEM identifier, because the medium of
   the main analysis treats vitamin B12 as a source it does not represent and stops if an exchange
   named EX_adocbl_e exists. A boundary reaction written in the import direction is reversed, so
   that a positive flux is export for every exchange, as in the main model.

2. Medium. src/medium_v2.py of the main analysis, unchanged, opens the same list of stems; a stem
   whose identifier Human-GEM's cross-reference does not provide is recorded as missing, as that
   module does for the main model. The biomass reaction of Human-GEM also requires cofactors that
   this list lacks. The inputs of Human-GEM's own growth task ("Growth on Ham's media",
   data/metabolicTasks/metabolicTasks_Essential.txt) that the list does not open are candidates; they
   are closed one at a time, in the order of the task file, and a candidate stays closed if biomass
   can still be produced. The candidates that have to stay open are added at an uptake capacity of
   1.0, the capacity of the assumed serum constituents and of Human-GEM's own growth check.

3. Reaction evidence. The gene-level matrix of the main analysis (data/gene_expression.npz, Entrez
   identifiers) is mapped to Human-GEM's genes through model/genes.tsv (a gene with several Entrez
   identifiers takes the larger value, the rule the main analysis uses for several probes of a gene),
   and the gene rules are evaluated by src/gpr.py with the conservative policy. The supported set is
   defined as in the main analysis: internal reactions with a rule and complete evidence in every
   retained profile.

4. configs/study.json names Human-GEM's biomass reaction as the objective. Every other setting is
   that of the main analysis.

The script writes second/exchange_map.tsv, second/medium_additions.tsv, second/build_report.json and
the root. It is deterministic.
"""
import argparse, csv, hashlib, json, shutil, sys
from pathlib import Path
import numpy as np
import cobra

BIOMASS = "MAR13082"
TASK_FILE = "data/metabolicTasks/metabolicTasks_Essential.txt"
TASK_ID = "GR"
ADDITION_CAPACITY = 1.0
GROWTH_FLOOR = 1e-6
COPY_DIRS = ("src", "mechanism", "closure")
COPY_FILES = ("environment.lock", "configs/study.json", "manifests/contexts.tsv", "manifests/metabolite_map.tsv",
              "manifests/analysis_panel.tsv", "protocol/analysis_plan.json", "protocol/PROTOCOL_LOCK.json",
              "data/gene_expression.npz")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


def task_inputs(path):
    """The IN metabolites (name[compartment]) of the growth task, as Human-GEM's own check reads them."""
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = rows[0]; id_col, in_col = header.index("ID"), header.index("IN")
    out, inside = [], False
    for row in rows[1:]:
        row = row + [""] * (len(header) - len(row))
        if row[id_col].strip():
            if inside: break
            inside = row[id_col].strip() == TASK_ID
        if inside:
            out += [t.strip() for t in row[in_col].split(";") if t.strip()]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--main", type=Path, required=True); ap.add_argument("--humangem", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True); ap.add_argument("--second", type=Path, default=Path(__file__).resolve().parent)
    a = ap.parse_args()
    if a.out.exists(): sys.exit(f"{a.out} exists; build into a new directory")
    hg = a.humangem
    rx_ann = {r["rxns"]: r for r in tsv(hg / "model/reactions.tsv")}
    met_ann = {r["mets"]: r for r in tsv(hg / "model/metabolites.tsv")}
    gene_ann = tsv(hg / "model/genes.tsv")
    model = cobra.io.read_sbml_model(str(hg / "model/Human-GEM.xml"))
    report = dict(humangem=dict(version=(hg / "version.txt").read_text().strip() if (hg / "version.txt").exists() else None,
                                files={f: sha(hg / f) for f in ("model/Human-GEM.xml", "model/reactions.tsv", "model/metabolites.tsv",
                                                                 "model/genes.tsv", TASK_FILE)}),
                  model=dict(reactions=len(model.reactions), metabolites=len(model.metabolites), genes=len(model.genes)))

    # ---- 1. boundary reactions renamed by Human-GEM's own cross-reference to BiGG
    boundary = [r for r in model.reactions if r.boundary]
    assert all(len(r.metabolites) == 1 for r in boundary)
    cand = {}
    for r in boundary:
        m = next(iter(r.metabolites)); big = (rx_ann.get(r.id, {}).get("rxnBiGGID") or "").strip()
        mb = (met_ann.get(m.id, {}).get("metBiGGID") or "").strip()
        if mb and m.compartment == "e": cand[r.id] = (f"EX_{mb}_e", "metBiGGID")
        elif big.startswith("EX_"): cand[r.id] = (big, "rxnBiGGID")
    count = {}
    for rid, (new, _) in cand.items(): count[new] = count.get(new, 0) + 1
    rows, reversed_n = [], 0
    for r in sorted(boundary, key=lambda x: x.id):
        m = next(iter(r.metabolites)); old = r.id; new, basis = cand.get(old, (None, None))
        if new is None or count[new] != 1: new, basis = f"EX_{old}", ("no cross-reference" if new is None else "cross-reference not unique")
        if new == "EX_adocbl_e": new, basis = f"EX_{old}", "kept: the main medium treats vitamin B12 as unrepresented"
        coef = r.metabolites[m]; flipped = False
        if coef > 0:                                   # written as import: reverse it, so that positive flux is export
            lb, ub = r.lower_bound, r.upper_bound
            r.add_metabolites({m: -2 * coef}); r.bounds = (-ub, -lb); flipped = True; reversed_n += 1
        rows.append(dict(humangem_reaction=old, pipeline_id=new, basis=basis, humangem_metabolite=m.id, metabolite_name=m.name,
                         compartment=m.compartment, metBiGGID=(met_ann.get(m.id, {}).get("metBiGGID") or ""),
                         metKEGGID=(met_ann.get(m.id, {}).get("metKEGGID") or ""), reversed=flipped,
                         lower_bound=r.lower_bound, upper_bound=r.upper_bound))
    new_ids = [x["pipeline_id"] for x in rows]
    assert len(set(new_ids)) == len(new_ids)
    for x in rows: model.reactions.get_by_id(x["humangem_reaction"]).id = x["pipeline_id"]
    model.repair()
    write_tsv(a.second / "exchange_map.tsv", rows, list(rows[0]))
    report["boundary"] = dict(n=len(rows), by_rxnBiGGID=sum(x["basis"] == "rxnBiGGID" for x in rows),
                              by_metBiGGID=sum(x["basis"] == "metBiGGID" for x in rows),
                              own_identifier=sum(x["pipeline_id"].startswith("EX_MAR") for x in rows), reversed=reversed_n,
                              compartments=sorted({x["compartment"] for x in rows}))

    # ---- targets and medium stems of the main analysis, as found in the renamed model
    mm = tsv(a.main / "manifests/metabolite_map.tsv"); plan = json.loads((a.main / "protocol/analysis_plan.json").read_text())
    by_new = {x["pipeline_id"]: x for x in rows}
    def kegg_ok(row, ex):
        k = {s.strip() for s in ex["metKEGGID"].replace(",", ";").split(";") if s.strip()}
        return bool(k) and (row.get("source_kegg", "") in k or row.get("model_kegg", "") in k)
    panel = [r for r in mm if str(r["primary_eligible"]).lower() in ("1", "true", "yes")]
    prim = {t["exchange_id"] for t in plan["primary_targets"]}
    tgt = []
    for r in panel:
        ex = by_new.get(r["exchange_id"])
        tgt.append(dict(metabolite_id=r["metabolite_id"], source_label=r["source_label"], exchange_id=r["exchange_id"],
                        primary=r["exchange_id"] in prim, present=ex is not None,
                        humangem_reaction=ex["humangem_reaction"] if ex else "", humangem_metabolite=ex["humangem_metabolite"] if ex else "",
                        metabolite_name=ex["metabolite_name"] if ex else "", basis=ex["basis"] if ex else "",
                        kegg_concordant=kegg_ok(r, ex) if ex else False, source_kegg=r.get("source_kegg", ""),
                        humangem_kegg=ex["metKEGGID"] if ex else "",
                        rxnBiGGID=(rx_ann.get(ex["humangem_reaction"], {}).get("rxnBiGGID") or "") if ex else "",
                        same_by_rxnBiGGID=bool(ex) and (rx_ann.get(ex["humangem_reaction"], {}).get("rxnBiGGID") or "") == r["exchange_id"]))
    write_tsv(a.second / "target_map.tsv", tgt, list(tgt[0]))
    report["targets"] = dict(panel=len(tgt), panel_present=sum(t["present"] for t in tgt), primary=sum(t["primary"] for t in tgt),
                             primary_present=sum(t["primary"] and t["present"] for t in tgt),
                             primary_absent=[t["exchange_id"] for t in tgt if t["primary"] and not t["present"]],
                             panel_absent=[t["exchange_id"] for t in tgt if not t["present"]],
                             primary_kegg_not_concordant=[t["exchange_id"] for t in tgt if t["primary"] and t["present"] and not t["kegg_concordant"]],
                             panel_not_same_by_rxnBiGGID=[t["exchange_id"] for t in tgt if t["present"] and not t["same_by_rxnBiGGID"]])

    # ---- the root: code and inputs of the main analysis, the model, the objective
    for d in COPY_DIRS: shutil.copytree(a.main / d, a.out / d, ignore=shutil.ignore_patterns("__pycache__"))
    for f in COPY_FILES:
        (a.out / f).parent.mkdir(parents=True, exist_ok=True); shutil.copy(a.main / f, a.out / f)
    cfg = json.loads((a.out / "configs/study.json").read_text()); cfg["objective"] = BIOMASS
    cfg["study_version"] = cfg["study_version"] + "+second-reconstruction"
    (a.out / "configs/study.json").write_text(json.dumps(cfg, indent=2) + "\n")
    (a.out / "data/raw").mkdir(parents=True, exist_ok=True)
    cobra.io.save_json_model(model, str(a.out / "data/raw/Recon3D.json"))

    # ---- 2. medium: the main list through the main module, then the growth-task inputs that have to stay open
    sys.path.insert(0, str(a.out))
    from src.medium_v2 import apply_medium
    def growth(extra):
        mdl = cobra.io.load_json_model(str(a.out / "data/raw/Recon3D.json")); mdl.solver = cfg["solver"]
        opened, missing, _ = apply_medium(mdl, {"assumed_serum": cfg["scenarios"]["primary"]["serum_uptake"]})
        for rid in extra: mdl.reactions.get_by_id(rid).lower_bound = -ADDITION_CAPACITY
        mdl.objective = BIOMASS; v = mdl.slim_optimize(error_value=float("nan"))
        return (0.0 if not np.isfinite(v) else float(v)), opened, missing
    g0, opened, missing = growth([])
    name_to_ex = {}
    for x in rows:
        if x["compartment"] == "e": name_to_ex.setdefault(x["metabolite_name"], []).append(x["pipeline_id"])
    cands, unresolved = [], []
    for token in task_inputs(hg / TASK_FILE):
        name, comp = token.rsplit("[", 1)
        ex = name_to_ex.get(name, []) if comp.rstrip("]") == "e" else []
        if len(ex) != 1: unresolved.append(token); continue
        if ex[0] not in opened and ex[0] not in [c[0] for c in cands]: cands.append((ex[0], name))
    g_all, _, _ = growth([c[0] for c in cands])
    keep = [c[0] for c in cands]; trail = []
    for rid, name in cands:
        trial = [k for k in keep if k != rid]; g, _, _ = growth(trial)
        stays_closed = g > GROWTH_FLOOR
        trail.append(dict(exchange_id=rid, metabolite_name=name, growth_without_it=g, added=not stays_closed))
        if stays_closed: keep = trial
    g_final, _, _ = growth(keep)
    add_rows = [dict(exchange_id=t["exchange_id"], metabolite_name=t["metabolite_name"], humangem_reaction=by_new[t["exchange_id"]]["humangem_reaction"],
                     added=t["added"], growth_without_it=repr(t["growth_without_it"]), uptake_capacity=ADDITION_CAPACITY if t["added"] else "",
                     basis="input of Human-GEM's growth task; biomass cannot be produced without it" if t["added"] else
                           "input of Human-GEM's growth task; not needed for biomass on this medium, left closed") for t in trail]
    write_tsv(a.second / "medium_additions.tsv", add_rows, list(add_rows[0]))
    report["medium"] = dict(opened_by_main_list=len(opened), missing_from_main_list=[m["exchange"] for m in missing],
                            growth_main_list_only=g0, task_inputs=len(task_inputs(hg / TASK_FILE)), task_inputs_unresolved=unresolved,
                            candidates=[c[0] for c in cands], growth_all_candidates=g_all, added=keep, growth_final=g_final)
    if not g_final > GROWTH_FLOOR: sys.exit("biomass cannot be produced on the final medium; see the report so far: " + json.dumps(report["medium"]))
    src = (a.out / "src/medium_v2.py").read_text()
    block = ('\n\n# ---- second reconstruction (written by second/build_root.py; the text above is the module of the main analysis, unchanged)\n'
             '# Inputs of Human-GEM\'s own growth task that the list above does not open and without which its biomass\n'
             '# reaction cannot carry flux. Uptake capacity as for the assumed serum constituents.\n'
             f'MODEL_GROWTH_TASK_INPUTS = {json.dumps([[rid, dict(cands)[rid]] for rid in keep])}\n'
             f'MODEL_GROWTH_TASK_CAPACITY = {ADDITION_CAPACITY!r}\n'
             '_apply_main_medium = apply_medium\n\n\n'
             'def apply_medium(model, uptake=None):\n'
             '    opened, missing, ledger = _apply_main_medium(model, uptake)\n'
             '    for rid, name in MODEL_GROWTH_TASK_INPUTS:\n'
             '        model.reactions.get_by_id(rid).lower_bound = -MODEL_GROWTH_TASK_CAPACITY\n'
             '        opened[rid] = {"lower_bound": -MODEL_GROWTH_TASK_CAPACITY, "category": "model_growth_task_input",\n'
             '                       "availability_basis": "input of the growth task of Human-GEM (Growth on Ham\'s media); needed by its biomass reaction",\n'
             '                       "bound_basis": "assumed transport capacity, as for the assumed serum constituents"}\n'
             '        ledger.append((rid, "lower_bound", 0.0, -MODEL_GROWTH_TASK_CAPACITY, "model_growth_task_input"))\n'
             '    return opened, missing, ledger\n')
    (a.out / "src/medium_v2.py").write_text(src + block)
    import importlib, src.medium_v2, src.cpu_model
    importlib.reload(src.medium_v2); importlib.reload(src.cpu_model)     # the pipeline now sees the extended medium
    from src.cpu_model import build_model

    # ---- 3. reaction evidence from the gene-level matrix of the main analysis
    from src.gpr import evaluate
    G = np.load(a.main / "data/gene_expression.npz", allow_pickle=False)
    genes = G["genes"].astype(str).tolist(); gi = {g: i for i, g in enumerate(genes)}; contexts = G["contexts"].astype(str).tolist()
    gvals, n_multi, gene_rows = {}, 0, []
    for row in gene_ann:
        ids = [s.strip() for s in (row.get("geneEntrezID") or "").replace(",", ";").split(";") if s.strip()]
        hit = [e for e in ids if e in gi]
        if hit:
            v = G["A"][gi[hit[0]]].copy()
            for e in hit[1:]: v = np.fmax(v, G["A"][gi[e]])
            gvals[row["genes"]] = v; n_multi += len(hit) > 1
        gene_rows.append(dict(model_gene_id=row["genes"], entrez_id=";".join(ids), measured=bool(hit)))
    write_tsv(a.out / "manifests/model_gene_map.tsv", gene_rows, list(gene_rows[0]))
    rx_ids, ev, cov = [], [], []
    for r in model.reactions:
        rule = (r.gene_reaction_rule or "").strip()
        if len(r.metabolites) <= 1 or not rule: continue
        v = evaluate(rule, gvals); ok = v is not None and np.isfinite(v).all()
        cov.append(dict(reaction_id=r.id, supported_all_expression_contexts=bool(ok)))
        if ok: rx_ids.append(r.id); ev.append(np.asarray(v, float))
    A = np.vstack(ev)
    np.savez_compressed(a.out / "data/scores_conservative.npz", A=A, rxn=np.array(rx_ids), contexts=np.array(contexts), cols=np.array(contexts))
    write_tsv(a.out / "manifests/gpr_coverage.tsv", cov, list(cov[0]))
    report["evidence"] = dict(model_genes=len(gene_ann), genes_measured=len(gvals), genes_with_several_measured_entrez=n_multi,
                              internal_reactions_with_rule=len(cov), supported_reactions=len(rx_ids), contexts=len(contexts),
                              evidence_min=float(A.min()), evidence_max=float(A.max()))

    # ---- the model as the pipeline builds it
    mdl, meta = build_model(a.out, cfg, "primary")
    report["pipeline_model"] = dict(g_max=meta["g_max"], task_flux=meta["task_flux"], objective=meta["objective"], objective_name=meta["objective_name"],
                                    opened=len(meta["opened"]), missing=[m["exchange"] for m in meta["missing"]])
    report["root_files"] = {f: sha(a.out / f) for f in ("data/raw/Recon3D.json", "data/scores_conservative.npz", "configs/study.json", "src/medium_v2.py",
                                                         "src/cpu_model.py", "src/cpu_readout.py", "src/run_study.py", "src/gpr.py")}
    report["unchanged_from_main"] = {f: sha(a.out / f) == sha(a.main / f) for f in ("src/cpu_model.py", "src/cpu_readout.py", "src/run_study.py", "src/gpr.py",
                                                                                    "src/evaluation.py", "manifests/contexts.tsv", "manifests/metabolite_map.tsv",
                                                                                    "protocol/analysis_plan.json")}
    (a.second / "build_report.json").write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps({k: report[k] for k in ("model", "boundary", "targets", "medium", "evidence", "pipeline_model")}, indent=1))


if __name__ == "__main__":
    main()
