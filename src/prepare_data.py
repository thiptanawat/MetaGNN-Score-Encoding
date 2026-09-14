"""Rebuild data from raw inputs, without outcome fitting or held-out evaluation.

Run from study root: .venv/bin/python -m src.prepare_data
Chemistry uses source/model annotation concordance, never observed CORE variation.
Only development observations determine replicate-noise tolerances.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import zipfile
from collections import Counter, defaultdict
import numpy as np
import xlrd
from .gpr import evaluate

DEFAULT_ROOT = Path(__file__).resolve().parents[1]
EXPRESSION_MEMBER = "output/RNA__Affy_HG_U133(A_B)_GCRMA.xls"
METADATA_MEMBER = "output/documentation/NCI60_CELL_LINE_METADATA.xls"
CHEMICAL_REFERENCE_SHA256 = "6512193df6635edba02ec46d4ecc758a633b2650048c060aec7dc552f199e752"
RAW_SOURCES = {
    "data/raw/Recon3D.json": ("aba925f17547a42f9fdb4c1f685d89364cbf4979bbe7862e9f793af7169b26d5", "http://bigg.ucsd.edu/static/models/Recon3D.json"),
    "data/raw/cellminer.zip": ("334dca10d425eaf9fa3c7fc386f8d413ef32d8408eb4fc5508ca75e6cf6b734e", "https://discover.nci.nih.gov/cellminer/download/processeddataset/nci60_RNA__Affy_HG_U133%28A_B%29_GCRMA.zip"),
    "data/raw/core/NIHMS419088-supplement-Database_S1.xls": ("9a089755b35e16fd17d61b50865ca695bbb33ec8bc8616b285cd5e1fdbfa2ce4", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3526189/"),
    "data/raw/core/NIHMS419088-supplement.pdf": ("1bf01aa2715af80924aa3dbc74c00507259f9517877f23e3157344966cc0f187", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3526189/"),
}
EXPOSED = {"RE:786-0", "LC:HOP-62", "LC:HOP-92", "BR:HS 578T", "CO:HT29", "ME:MALME-3M", "BR:MDA-MB-231", "LC:NCI-H226", "LE:RPMI-8226", "LE:SR", "RE:UO-31"}
ORIGINS = {**dict.fromkeys(["ME:M14", "ME:MDA-MB-435", "ME:MDA-N"], "M14_MDA-MB-435"),
           **dict.fromkeys(["OV:OVCAR-8", "OV:NCI/ADR-RES"], "OVCAR-8_ADR-RES"),
           **dict.fromkeys(["CNS:SNB-19", "CNS:U251"], "SNB-19_U251")}
LINEAGES = {"BR": "breast", "CNS": "central nervous system", "CO": "colon", "LC": "lung", "LE": "leukemia", "ME": "melanoma", "OV": "ovary", "PR": "prostate", "RE": "renal"}
SOURCE_ALIASES = {
    "DHAP": "dihydroxyacetone phosphate", "2-aminodipate": "2-aminoadipate",
    "PEP": "phosphoenolpyruvate", "ADMA": "asymmetric dimethylarginine", "NMMA": "N-monomethylarginine",
    "fru-1,6-DP/fru-2,6-DP/glc-1,6-DP": "fructose-1,6-diphosphate/fructose-2,6-diphosphate/glucose-1,6-diphosphate",
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(rows[0])
    path = Path(path)
    staging = path.with_name(path.name + ".preparation-tmp")
    with staging.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    staging.replace(path)


def write_json(path, obj):
    path = Path(path)
    staging = path.with_name(path.name + ".preparation-tmp")
    staging.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n")
    staging.replace(path)


def normalize_context(label):
    label = str(label).strip().upper()
    label = re.sub(r"/ATCC$|/TB$|\(TB\)$", "", label)
    label = re.sub(r"[^A-Z0-9]", "", label)
    return "7860" if label == "786O" else label


def unique_normalized(labels, expression=False):
    result = {}
    for label in labels:
        key = normalize_context(label.split(":", 1)[1] if expression else label)
        if key in result:
            raise ValueError(f"Ambiguous context normalization: {label}, {result[key]}")
        result[key] = label
    return result


def finite_number(value):
    try:
        number = float(value)
    except (ValueError, TypeError):
        return np.nan
    return number if np.isfinite(number) else np.nan


def parse_expression(root, model):
    with zipfile.ZipFile(root / "data/raw/cellminer.zip") as archive:
        contents = archive.read(EXPRESSION_MEMBER)
        metadata_contents = archive.read(METADATA_MEMBER)
        members = [{"member": name, "bytes": len(archive.read(name)), "sha256": hashlib.sha256(archive.read(name)).hexdigest()} for name in sorted(archive.namelist())]
    write_tsv(root / "manifests/expression_archive_members.tsv", members)
    sheet = xlrd.open_workbook(file_contents=contents).sheet_by_index(0)
    header_row = next(row for row in range(20) if str(sheet.cell_value(row, 0)).startswith("Identifier"))
    headers = [str(value).strip() for value in sheet.row_values(header_row)]
    assert headers[2].startswith("Entrez gene id")
    columns = headers[7:]
    probe_rows, probe_values, valid_entrez = [], [], []
    for row in range(header_row + 1, sheet.nrows):
        values = np.array([finite_number(sheet.cell_value(row, col)) for col in range(7, sheet.ncols)])
        gene_number = finite_number(sheet.cell_value(row, 2))
        gene = str(int(gene_number)) if np.isfinite(gene_number) and gene_number > 0 and gene_number.is_integer() else ""
        probe_rows.append({"source_excel_row": row + 1, "probe_id": str(sheet.cell_value(row, 0)), "entrez_id": gene,
                           "gene_symbol": str(sheet.cell_value(row, 1)), "included": bool(gene), "reason": "positive Entrez identifier" if gene else "no positive Entrez identifier"})
        if gene:
            valid_entrez.append(gene)
            probe_values.append(values)
    matrix = np.vstack(probe_values)
    keep = np.isfinite(matrix).any(axis=0)
    retained = [column for column, flag in zip(columns, keep) if flag]
    matrix = matrix[:, keep]
    genes = {}
    for gene, values in zip(valid_entrez, matrix):
        genes[gene] = values.copy() if gene not in genes else np.fmax(genes[gene], values)
    gene_ids = sorted(genes, key=int)
    gene_matrix = np.vstack([genes[gene] for gene in gene_ids])
    model_genes = {}
    for gene in model["genes"]:
        base = gene["id"].split(".")[0].split("_")[0]
        model_genes[gene["id"]] = base if base.isdigit() and base != "0" else ""
    gvals = {identifier: genes[base] for identifier, base in model_genes.items() if base in genes}
    coverage, reaction_ids, evidence = [], [], []
    for reaction in model["reactions"]:
        rule = reaction.get("gene_reaction_rule", "").strip()
        if len(reaction["metabolites"]) <= 1 or not rule:
            continue
        value = evaluate(rule, gvals)
        supported = value is not None and np.isfinite(value).all()
        coverage.append({"reaction_id": reaction["id"], "gpr": rule, "supported_all_expression_contexts": bool(supported),
                         "missing_context_count": len(retained) if value is None else int((~np.isfinite(value)).sum()),
                         "reason": "complete measured GPR alternative" if supported else "no complete measured alternative in every retained context"})
        if supported:
            reaction_ids.append(reaction["id"])
            evidence.append(value)
    scores = np.vstack(evidence)
    np.savez_compressed(root / "data/scores_conservative.npz", A=scores, rxn=np.array(reaction_ids), contexts=np.array(retained), cols=np.array(retained))
    np.savez_compressed(root / "data/gene_expression.npz", A=gene_matrix, genes=np.array(gene_ids), contexts=np.array(retained))
    write_tsv(root / "manifests/probe_gene_map.tsv", probe_rows)
    write_tsv(root / "manifests/gpr_coverage.tsv", coverage)
    write_tsv(root / "manifests/model_gene_map.tsv", [{"model_gene_id": identifier, "entrez_id": base, "measured": base in genes} for identifier, base in model_genes.items()])
    metadata = xlrd.open_workbook(file_contents=metadata_contents).sheet_by_index(0)
    tissue = {}
    for row in range(8, metadata.nrows):
        value = str(metadata.cell_value(row, 0)).strip()
        if ":" in value:
            tissue[normalize_context(value.split(":", 1)[1])] = str(metadata.cell_value(row, 1)).strip()
    summary = {"archive_member": EXPRESSION_MEMBER, "header_excel_row": header_row + 1, "database_version": str(sheet.cell_value(2, 1)),
               "genome_version": str(sheet.cell_value(3, 1)), "source_date_literal": str(sheet.cell_value(4, 1)),
               "source_scale": "CellMiner GCRMA intensity, not TPM or measured protein", "aggregation": "per-context maximum of probes with the same positive Entrez ID",
               "raw_probe_rows": len(probe_rows), "positive_entrez_probe_rows": len(valid_entrez), "unique_entrez": len(genes),
               "negative_probe_values_retained_without_clipping": int((matrix < 0).sum()),
               "expression_columns": len(columns), "retained_columns": len(retained), "dropped_all_missing_columns": [column for column, flag in zip(columns, keep) if not flag],
               "retained_gene_missing_cells": int((~np.isfinite(gene_matrix)).sum()), "model_base_entrez": len(set(model_genes.values()) - {""}),
               "measured_model_base_entrez": len(set(model_genes.values()) & set(genes)), "supported_reactions": len(reaction_ids),
               "gpr_policy": "AND minimum propagating unknown; OR maximum over supported alternatives; fixed fully finite reaction set"}
    write_json(root / "data/expression_summary.json", summary)
    return columns, retained, tissue, summary


def parse_table_s1(root):
    pdf = root / "data/raw/core/NIHMS419088-supplement.pdf"
    output = root / "data/source_supplement_layout.txt"
    subprocess.run(["pdftotext", "-layout", str(pdf), str(output)], check=True)
    text = output.read_text()
    records, active = [], False
    for page_number, page in enumerate(text.split("\f"), 1):
        if "Table S1. Monitored metabolites" in page:
            active = True
        if not active:
            continue
        previous = ""
        for line in page.splitlines():
            method = re.search(r"\b(HILIC|IPR|BGA)\b", line)
            kegg = re.search(r"((?:C\d{5}|NA)(?:/(?:C\d{5}|NA))*)\s*(?:\s+[01]){0,2}\s*$", line)
            if method and kegg and kegg.start() > method.end():
                name = line[:method.start()].strip() or previous
                if name and "Method" not in name:
                    records.append({"source_label": name, "method": method.group(), "source_kegg": kegg.group(1), "source_pdf_page": page_number,
                                    "source_text": line.strip() if line[:method.start()].strip() else previous + " " + line.strip()})
            if line.strip():
                previous = line.strip()
    if len(records) != 219 or len({(row["source_label"], row["method"]) for row in records}) != 219:
        raise ValueError("Table S1 extraction contract failed: expected 219 unique method/analyte records")
    write_tsv(root / "manifests/source_table_s1.tsv", records)
    return {(row["source_label"], row["method"]): row for row in records}


def parse_core(root):
    book = xlrd.open_workbook(root / "data/raw/core/NIHMS419088-supplement-Database_S1.xls")
    core = book.sheet_by_name("CORE data")
    groups = defaultdict(list)
    for column in range(3, core.ncols):
        groups[str(core.cell_value(0, column)).strip()].append(column)
    if len(groups) != 60 or {len(cols) for cols in groups.values()} != {2}:
        raise ValueError("CORE replicate structure changed")
    peaks = book.sheet_by_name("MS peak areas")
    peak_headers = [str(value).strip() for value in peaks.row_values(0)]
    historical = {}
    for row in range(1, peaks.nrows):
        historical[str(peaks.cell_value(row, 1)).strip()] = {
            "historical_variable": bool(peaks.cell_value(row, peak_headers.index("Variable (d)")) == 1),
            "historical_sd_ratio": finite_number(peaks.cell_value(row, peak_headers.index("SD Ratio (c)")))}
    metabolites, values = [], {}
    for row in range(1, core.nrows):
        name = str(core.cell_value(row, 1)).strip()
        metabolite_id = f"core_{row:03d}"
        metabolites.append({"metabolite_id": metabolite_id, "source_excel_row": row + 1, "source_label": name,
                            "method": str(core.cell_value(row, 0)).strip(), "calibrated": bool(core.cell_value(row, 2) == 1), **historical[name]})
        values[metabolite_id] = {label: [finite_number(core.cell_value(row, col)) for col in cols] for label, cols in groups.items()}
    if len(metabolites) != 140 or not all(np.isfinite(v).all() for group in values.values() for v in group.values()):
        raise ValueError("CORE source shape or numeric completeness changed")
    return metabolites, values, list(groups)


def make_contexts(root, raw_columns, retained, tissues, core_labels):
    expr = unique_normalized(raw_columns, expression=True)
    core = unique_normalized(core_labels)
    rows, core_to_context = [], {}
    for key in sorted(set(expr) | set(core)):
        expression, core_label = expr.get(key, ""), core.get(key, "")
        context = expression or ("BR:MDA-MB-468" if key == "MDAMB468" else "CORE:" + core_label)
        present_expression = expression in retained
        matched = present_expression and bool(core_label)
        partition = "development" if matched and context in EXPOSED else "test" if matched else "excluded"
        reason = "previous source-label exposure; development only" if partition == "development" else "matched profile, not in 11-line exposure list" if matched else "no CORE outcome" if not core_label else "all-missing expression" if expression else "no expression profile"
        rows.append({"context_id": context, "origin_group": ORIGINS.get(context, context), "partition": partition,
                     "lineage": LINEAGES.get(context.split(":", 1)[0], tissues.get(key, "unknown").lower()),
                     "source_metadata_tissue": tissues.get(key, ""),
                     "normalized_id": key, "expression_column": expression, "core_label": core_label,
                     "has_expression": present_expression, "has_core": bool(core_label), "reason": reason})
        if core_label:
            core_to_context[core_label] = context
    development = [row["context_id"] for row in rows if row["partition"] == "development"]
    test = [row["context_id"] for row in rows if row["partition"] == "test"]
    origins = {row["origin_group"] for row in rows if row["partition"] == "test"}
    if set(development) != EXPOSED or len(test) != 47 or len(origins) != 44:
        raise ValueError("Context identity contract failed")
    write_tsv(root / "manifests/contexts.tsv", rows)
    return rows, core_to_context


def make_chemistry_map(root, metabolites, table, model):
    model_mets = {met["id"]: met for met in model["metabolites"]}
    kegg_exchanges = defaultdict(list)
    for reaction in model["reactions"]:
        if not reaction["id"].startswith("EX_") or len(reaction["metabolites"]) != 1:
            continue
        met_id, coefficient = next(iter(reaction["metabolites"].items()))
        met = model_mets[met_id]
        if met.get("compartment") != "e" or coefficient == 0:
            continue
        kegg_ids = met.get("annotation", {}).get("kegg.compound", [])
        if isinstance(kegg_ids, str):
            kegg_ids = [kegg_ids]
        for kegg in kegg_ids:
            kegg_exchanges[kegg].append((reaction, met, coefficient))
    rows = []
    for item in metabolites:
        name, method = item["source_label"], item["method"]
        source_name = SOURCE_ALIASES.get(name, name)
        duplicate_glycerol = name in {"glycerol_1", "glycerol_2"}
        source = table.get(("glycerol" if duplicate_glycerol else source_name, method))
        if source is None:
            raise ValueError(f"Unresolved source Table S1 label: {name}, {method}")
        kegg = source["source_kegg"]
        candidates = kegg_exchanges.get(kegg, [])
        row = {**item, "table_s1_label": "glycerol / glycerol_early (assay correspondence unresolved)" if duplicate_glycerol else source_name,
               "source_pdf_page": source["source_pdf_page"], "source_kegg": kegg,
               "exchange_id": "", "model_metabolite_id": "", "model_metabolite_name": "", "model_kegg": "", "model_formula": "", "sign_factor": "",
               "candidate_exchanges": ";".join(candidate[0]["id"] for candidate in candidates),
               "status": "excluded", "primary_eligible": False, "reason": ""}
        if "/" in name or "/" in kegg:
            row["reason"] = "source assay aggregates unresolved analytes"
        elif duplicate_glycerol:
            row["reason"] = "two glycerol assay labels share one chemical identifier; source assay correspondence unresolved; both excluded"
        elif name == "homocystine":
            row["reason"] = "source label/KEGG conflict: homocystine disulfide labeled C00155, the homocysteine identifier"
        elif name == "3-OH-kynurenate":
            row["reason"] = "source label/KEGG conflict: 3-OH-kynurenate labeled C03227, the 3-hydroxy-L-kynurenine identifier"
        elif name == "creatinine":
            row["reason"] = "model annotation conflict: creat_e has creatinine KEGG C00791 but creatine name/formula C4H9N3O2; no ad hoc remapping"
        elif kegg == "NA":
            row["reason"] = "no source KEGG identifier for auditable exact mapping"
        elif len(candidates) == 0:
            row["reason"] = "no extracellular EX reaction with the exact source KEGG identifier"
        elif len(candidates) > 1:
            row["reason"] = "multiple extracellular EX reactions share source KEGG identifier; no arbitrary selection"
        else:
            reaction, met, coefficient = candidates[0]
            row.update(exchange_id=reaction["id"], model_metabolite_id=met["id"], model_metabolite_name=met["name"],
                       model_kegg=";".join(met.get("annotation", {}).get("kegg.compound", [])), model_formula=met.get("formula", ""),
                       sign_factor=float(-coefficient), status="mapped", primary_eligible=item["calibrated"],
                       reason="single source/model KEGG-concordant extracellular exchange" if item["calibrated"] else "chemical mapping retained; uncalibrated assay excluded from primary QC panel")
        rows.append(row)
    counts = Counter(row["exchange_id"] for row in rows if row["primary_eligible"])
    for row in rows:
        if row["primary_eligible"] and counts[row["exchange_id"]] > 1:
            row.update(status="excluded", primary_eligible=False, reason="multiple eligible assays target the same exchange; all excluded pending external assay clarification")
    write_tsv(root / "manifests/metabolite_map.tsv", rows)
    return rows


def make_outcomes_and_noise(root, metabolites, values, context_rows, core_to_context):
    outcomes = []
    for item in metabolites:
        unit = "fmol/cell/h" if item["calibrated"] else "arbitrary"
        for label, pair in values[item["metabolite_id"]].items():
            for replicate, value in enumerate(pair, 1):
                outcomes.append({"context_id": core_to_context[label], "metabolite_id": item["metabolite_id"], "replicate": replicate, "value": repr(value), "unit": unit,
                                 "core_label": label, "source_excel_row": item["source_excel_row"], "calibrated": item["calibrated"], "historical_variable": item["historical_variable"]})
    write_tsv(root / "data/outcomes_long.tsv", outcomes)
    dev_labels = [row["core_label"] for row in context_rows if row["partition"] == "development"]
    noise = []
    for item in metabolites:
        differences = [abs(values[item["metabolite_id"]][label][0] - values[item["metabolite_id"]][label][1]) for label in dev_labels]
        noise.append({"metabolite_id": item["metabolite_id"], "tolerance": float(np.quantile(differences, .90, method="linear")),
                      "n_dev_pairs": len(differences), "rule": "90th percentile absolute replicate difference; 11 development contexts; numpy linear quantile", "unit": "fmol/cell/h" if item["calibrated"] else "arbitrary"})
    write_tsv(root / "data/assay_noise_development.tsv", noise)
    return len(outcomes)


def verify_chemical_reference(root, chemistry):
    """Check frozen primary-database facts without requiring a live service.

Hydrogen differences can arise from model protonation and are not treated as a
compound mismatch. Matching nonhydrogen formulas is necessary, not sufficient,
for identity: source/model identifiers and the explicit name-conflict exclusions
remain part of the mapping rule. The reference fixture contains facts, not full
KEGG records, and records primary URL, access date and captured-response hash.
"""
    path = root / "protocol/chemical_identity_reference.tsv"
    if digest(path) != CHEMICAL_REFERENCE_SHA256:
        raise ValueError("Frozen chemical reference checksum mismatch")
    with path.open() as handle:
        reference = {row["kegg_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    def nonhydrogen(formula):
        tokens = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
        if "".join(element + number for element, number in tokens) != formula:
            raise ValueError(f"Unsupported chemical formula: {formula}")
        return {element: int(number or 1) for element, number in tokens if element != "H"}
    review = []
    for row in chemistry:
        if row["status"] != "mapped":
            continue
        identity = reference[row["source_kegg"]]
        concordant = nonhydrogen(identity["formula"]) == nonhydrogen(row["model_formula"])
        if not concordant:
            raise ValueError(f"Accepted chemical mapping has elemental conflict: {row['metabolite_id']}")
        review.append({"metabolite_id": row["metabolite_id"], "source_label": row["source_label"], "source_kegg": row["source_kegg"],
                       "kegg_name": identity["name"], "kegg_formula": identity["formula"], "model_formula": row["model_formula"],
                       "nonhydrogen_formula_match": concordant, "reference_url": identity["source_url"], "access_date": identity["access_date"],
                       "primary_eligible": row["primary_eligible"]})
    write_tsv(root / "manifests/chemical_identity_review.tsv", review)
    write_json(root / "manifests/chemical_reference_provenance.json", {"path": str(path.relative_to(root)), "sha256": CHEMICAL_REFERENCE_SHA256,
               "entries": len(reference), "mapped_assays_checked": len(review), "source": "KEGG COMPOUND entries; primary database factual curation accessed 2026-09-13",
               "scope": "Identifier/name/formula audit; not independent structural identification or proof of chromatographic separation"})


def main(root=DEFAULT_ROOT):
    root = Path(root).resolve()
    for folder in ["data", "manifests", "reports"]:
        (root / folder).mkdir(parents=True, exist_ok=True)
    sources = []
    for relative, (expected, url) in RAW_SOURCES.items():
        actual = digest(root / relative)
        if actual != expected:
            raise ValueError(f"Raw source SHA256 mismatch: {relative}")
        sources.append({"path": relative, "sha256": actual, "bytes": (root / relative).stat().st_size, "source_url": url})
    write_tsv(root / "manifests/raw_sources.tsv", sources)
    model = json.loads((root / "data/raw/Recon3D.json").read_text())
    raw_columns, retained, tissues, expression = parse_expression(root, model)
    source_table = parse_table_s1(root)
    metabolites, values, core_labels = parse_core(root)
    context_rows, core_to_context = make_contexts(root, raw_columns, retained, tissues, core_labels)
    chemistry = make_chemistry_map(root, metabolites, source_table, model)
    verify_chemical_reference(root, chemistry)
    outcome_count = make_outcomes_and_noise(root, metabolites, values, context_rows, core_to_context)
    summary = {"schema_version": 1, "sources": sources, "expression": expression,
               "raw_core_contexts": len(core_labels), "outcome_rows": outcome_count, "source_table_s1_records": len(source_table),
               "development_contexts": [row["context_id"] for row in context_rows if row["partition"] == "development"],
               "test_contexts": [row["context_id"] for row in context_rows if row["partition"] == "test"],
               "test_origin_groups": sorted({row["origin_group"] for row in context_rows if row["partition"] == "test"}),
               "primary_chemistry_metabolite_ids": [row["metabolite_id"] for row in chemistry if row["primary_eligible"]],
               "primary_chemistry_count": sum(row["primary_eligible"] for row in chemistry),
               "mapping_status_counts": dict(Counter(row["status"] for row in chemistry)),
               "historical_variable_used_for_primary": False,
               "primary_panel_rule": "calibrated; single source chemical identifier; one extracellular model exchange; no ambiguous assays, identity conflict or duplicate target; later development-only observed estimability",
               "sign_convention": "CORE positive release; multiply model exchange flux by -stoichiometric coefficient of extracellular metabolite",
               "outcome_scale_caveat": "Absolute model flux units are not converted to fmol/cell/h; only predeclared within-metabolite comparisons may be evaluated",
               "noise_rule": "per-metabolite 90th percentile absolute replicate difference in 11 development contexts; numpy linear quantile",
               "no_test_performance_evaluation": True}
    write_json(root / "data/preparation_summary.json", summary)
    print(json.dumps({"reaction_evidence_shape": [expression["supported_reactions"], expression["retained_columns"]], "development": len(summary["development_contexts"]), "test_profiles": len(summary["test_contexts"]), "test_origins": len(summary["test_origin_groups"]), "chemistry_panel": summary["primary_chemistry_count"], "raw_outcome_rows": outcome_count}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    main(parser.parse_args().root)
