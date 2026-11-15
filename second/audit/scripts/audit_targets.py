"""Audit the 96 measured targets: main-model metabolite vs Human-GEM metabolite vs measured analyte.

Usage: python -I audit_targets.py OUTDIR
Writes target_audit.tsv and alternatives.tsv (long format) into OUTDIR.
"""
import csv, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_lib import *

out = Path(sys.argv[1])
M = Models()
targets = read_tsv(SEC / "target_map.tsv")
mapfile = {r["metabolite_id"]: r for r in read_tsv(MAPFILE)}

hg_ex, hg_e_info = build_hg_ex(M)

rows, alt_rows = [], []
for t in targets:
    mf = mapfile[t["metabolite_id"]]
    exid = t["exchange_id"]
    sec_rx = M.sr.get(exid); main_rx = M.mr.get(exid)
    row = dict(metabolite_id=t["metabolite_id"], source_label=t["source_label"], exchange_id=exid, primary=t["primary"],
               source_kegg=t["source_kegg"])
    main_met = mf["model_metabolite_id"]
    mi = M.main_info(main_met)
    hid = t["humangem_metabolite"]
    hi = M.hg_info(hid)
    # exchange checks
    if sec_rx:
        mets = sec_rx["metabolites"]
        row.update(ex_n_metabolites=len(mets), ex_metabolite=";".join(mets), ex_compartment=";".join(M.sm[m]["compartment"] for m in mets),
                   ex_coefficient=";".join(str(v) for v in mets.values()), ex_lower=sec_rx["lower_bound"], ex_upper=sec_rx["upper_bound"],
                   ex_matches_target_map=str(list(mets) == [hid]),
                   ex_reaction_in_map=M.ex_by_id.get(exid, {}).get("humangem_reaction", ""),
                   ex_rule_ok=str(len(mets) == 1 and all(M.sm[m]["compartment"] == "e" for m in mets) and list(mets.values()) == [-1.0]))
    else:
        row.update(ex_rule_ok="MISSING")
    if main_rx:
        mm_ = main_rx["metabolites"]
        row.update(main_ex_metabolite=";".join(mm_), main_ex_coefficient=";".join(str(v) for v in mm_.values()),
                   main_ex_lower=main_rx["lower_bound"], main_ex_upper=main_rx["upper_bound"],
                   main_ex_metabolite_matches_map=str(list(mm_) == [main_met]))
    else:
        row.update(main_ex_metabolite="MISSING")
    row.update(main_metabolite=main_met, hg_reaction=t["humangem_reaction"], hg_metabolite=hid, mapping_basis=t["basis"],
               main_name=mi["name"] if mi else "", hg_name=hi["name"] if hi else "",
               main_formula=mi["formula"] if mi else "", hg_formula=hi["formula"] if hi else "",
               main_formula_in_map=mf["model_formula"])
    cls, d = classify(mi, hi) if (mi and hi) else ("CANNOT_TELL", {"missing": "metabolite absent"})
    row.update(d)
    sk = set(split_ids(t["source_kegg"])); mk = set(mi["kegg"]) if mi else set(); hk = set(hi["kegg"]) if hi else set()
    row.update(main_kegg=";".join(mi["kegg"]) if mi else "", hg_kegg=";".join(hi["kegg"]) if hi else "",
               kegg_src_in_main=str(bool(sk & mk)), kegg_src_in_hg=str(bool(sk & hk)),
               kegg_concordant_target_map=t["kegg_concordant"], rxnBiGGID_same=t["same_by_rxnBiGGID"], rxnBiGGID=t["rxnBiGGID"],
               main_chebi=";".join(mi["chebi"]) if mi else "", hg_chebi=";".join(hi["chebi"]) if hi else "",
               main_hmdb=";".join(mi["hmdb"]) if mi else "", hg_hmdb=";".join(hi["hmdb"]) if hi else "",
               hg_inchi=hi.get("inchi", "") if hi else "", hg_smiles=hi.get("smiles", "") if hi else "",
               hg_bigg_ids=";".join(hi["bigg"]) if hi else "")
    row["class"] = cls
    al = alternatives(M, hg_ex, hg_e_info, hid, mi, t["source_kegg"]) if hid in hg_e_info else []
    strong = [a for a in al if a[2]]
    kshared = set(split_ids(t["source_kegg"])) | (set(mi["kegg"]) if mi else set())
    noex = [m for m, t_ in M.hgt.items() if m.endswith("e") and M.sm.get(m, {}).get("compartment") == "e" and m not in hg_ex
            and set(split_ids(t_["metKEGGID"])) & kshared]
    row["e_metabolites_without_exchange_sharing_kegg"] = ";".join("%s(%s)" % (m, M.sm[m]["name"]) for m in noex)
    row["n_alt_strong"] = len(strong); row["n_alt_weak"] = len(al) - len(strong)
    row["alt_strong"] = ";".join("%s(%s|%s)" % (a[0], ",".join(hg_ex[a[0]]), a[1]["name"]) for a in strong)
    for a in al:
        alt_rows.append(dict(metabolite_id=t["metabolite_id"], source_label=t["source_label"], mapped_exchange=exid,
                             mapped_hg_metabolite=hid, mapped_hg_name=hi["name"], alt_hg_metabolite=a[0], alt_name=a[1]["name"],
                             alt_formula=a[1]["formula"], alt_charge=a[1]["charge"], alt_exchange=",".join(hg_ex[a[0]]),
                             alt_exchange_bounds=";".join("%s,%s" % (M.sr[e]["lower_bound"], M.sr[e]["upper_bound"]) for e in hg_ex[a[0]]),
                             alt_kegg=";".join(a[1]["kegg"]), alt_chebi=";".join(a[1]["chebi"]), alt_inchikey=";".join(a[1]["keys"]),
                             strong_tags=a[2], weak_tags=a[3]))
    rows.append(row)

def wr(path, rs):
    fields = []
    for r in rs:
        for k in r:
            if k not in fields: fields.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader(); w.writerows(rs)

wr(out / "target_audit.tsv", rows)
wr(out / "alternatives.tsv", alt_rows)
import collections
print(collections.Counter(r["class"] for r in rows))
print("alt rows", len(alt_rows), "strong rows", sum(1 for a in alt_rows if a["strong_tags"]))
