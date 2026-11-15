"""Audit the medium inputs: every stem of BASAL, CULTURE and SERUM in src/medium_v2.py, and the six inputs
added by the second-root build (second/medium_additions.tsv).

Usage: python -I audit_medium.py OUTDIR
Writes medium_audit.tsv, medium_missing_candidates.tsv and medium_alternatives.tsv into OUTDIR.
"""
import csv, re, sys, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_lib import *

out = Path(sys.argv[1])
M = Models()
hg_ex, hg_e_info = build_hg_ex(M)
ML = medium_lists()
adds = read_tsv(SEC / "medium_additions.tsv")
added_rows = [a for a in adds if a["added"] == "True"]
opened_main_stems = set()

# main-model extracellular metabolites that have a boundary exchange
main_ex = {}
for r in M.main["reactions"]:
    if len(r["metabolites"]) == 1:
        mid = next(iter(r["metabolites"]))
        if M.mm[mid]["compartment"] == "e":
            main_ex.setdefault(mid, []).append(r["id"])


def main_counterpart(stem, hg_inf):
    """Main-model metabolite for a stem: by its exchange, else by the id stem_e, else by shared KEGG or InChIKey skeleton."""
    rid = "EX_%s_e" % stem
    if rid in M.mr:
        return next(iter(M.mr[rid]["metabolites"])), "main exchange EX_%s_e" % stem
    if "%s_e" % stem in M.mm:
        return "%s_e" % stem, "main metabolite id %s_e (no exchange)" % stem
    if hg_inf:
        for mid, m in M.mm.items():
            if m["compartment"] != "e": continue
            inf = M.main_info(mid)
            if set(inf["kegg"]) & set(hg_inf["kegg"]):
                return mid, "shared KEGG"
    return None, ""


def search_main(hg_inf):
    """Main-model extracellular metabolites with an exchange that could be the same compound as hg_inf."""
    res = []
    hk = {key_parts(k)[0] for k in hg_inf["keys"] if key_parts(k)}
    for mid, exs in main_ex.items():
        inf = M.main_info(mid)
        tags = []
        if set(inf["kegg"]) & set(hg_inf["kegg"]): tags.append("kegg:" + ",".join(sorted(set(inf["kegg"]) & set(hg_inf["kegg"]))))
        if set(inf["chebi"]) & set(hg_inf["chebi"]): tags.append("chebi")
        if {key_parts(k)[0] for k in inf["keys"] if key_parts(k)} & hk: tags.append("inchikey_skeleton")
        if norm_name(inf["name"]) == norm_name(hg_inf["name"]): tags.append("name")
        if hg_inf["bigg"] and set(inf["bigg"]) & set(hg_inf["bigg"]): tags.append("bigg")
        if tags:
            res.append((mid, exs, inf, ";".join(tags)))
    return res

medium_stems = []
for lst in ("BASAL", "CULTURE", "SERUM"):
    for stem, cat in ML[lst]:
        medium_stems.append((stem, cat, lst))

rows, alt_rows, cand_rows = [], [], []


def add_ex_checks(row, rid, hg_met_expected=None):
    r = M.sr.get(rid)
    if not r:
        row["ex_rule_ok"] = "NOT_IN_SECOND_ROOT"; return None
    mets = r["metabolites"]
    row.update(ex_n_metabolites=len(mets), ex_metabolite=";".join(mets), ex_compartment=";".join(M.sm[m]["compartment"] for m in mets),
               ex_coefficient=";".join(str(v) for v in mets.values()), ex_lower=r["lower_bound"], ex_upper=r["upper_bound"],
               ex_rule_ok=str(len(mets) == 1 and all(M.sm[m]["compartment"] == "e" for m in mets) and list(mets.values()) == [-1.0]),
               ex_reaction_in_map=M.ex_by_id.get(rid, {}).get("humangem_reaction", ""),
               ex_basis=M.ex_by_id.get(rid, {}).get("basis", ""))
    return next(iter(mets)) if len(mets) == 1 else None


def fill(row, mi, hi):
    cls, d = classify(mi, hi) if (mi and hi) else ("CANNOT_TELL", {"missing": "metabolite absent"})
    row.update(d)
    row.update(main_name=mi["name"] if mi else "", hg_name=hi["name"] if hi else "",
               main_formula=mi["formula"] if mi else "", hg_formula=hi["formula"] if hi else "",
               main_kegg=";".join(mi["kegg"]) if mi else "", hg_kegg=";".join(hi["kegg"]) if hi else "",
               main_chebi=";".join(mi["chebi"]) if mi else "", hg_chebi=";".join(hi["chebi"]) if hi else "",
               main_hmdb=";".join(mi["hmdb"]) if mi else "", hg_hmdb=";".join(hi["hmdb"]) if hi else "",
               hg_inchi=hi.get("inchi", "") if hi else "", hg_bigg_ids=";".join(hi["bigg"]) if hi else "")
    row["class"] = cls


for stem, cat, lst in medium_stems:
    rid = "EX_%s_e" % stem
    row = dict(kind="medium_stem", list=lst, stem=stem, category=cat, exchange_id=rid,
               in_main_model=str(rid in M.mr), in_second_root=str(rid in M.sr))
    row["opened_by_main_list"] = str(rid in M.mr and cat != "unrepresented_source_b12")
    if rid in M.mr:
        mr = M.mr[rid]
        row.update(main_ex_metabolite=";".join(mr["metabolites"]), main_ex_coefficient=";".join(str(v) for v in mr["metabolites"].values()),
                   main_ex_lower=mr["lower_bound"], main_ex_upper=mr["upper_bound"])
    hid = add_ex_checks(row, rid)
    row["hg_metabolite"] = hid or ""
    hi = M.hg_info(hid) if hid else None
    mid, how = main_counterpart(stem, hi)
    row["main_metabolite"] = mid or ""; row["main_counterpart_basis"] = how
    mi = M.main_info(mid) if mid else None
    if hi:
        fill(row, mi, hi)
        if not mi: row["missing"] = "no main-model metabolite found for this stem"
        al = alternatives(M, hg_ex, hg_e_info, hid, mi, "")
        strong = [a for a in al if a[2]]
        row["n_alt_strong"] = len(strong); row["n_alt_weak"] = len(al) - len(strong)
        row["alt_strong"] = ";".join("%s(%s|%s)" % (a[0], ",".join(hg_ex[a[0]]), a[1]["name"]) for a in strong)
        for a in al:
            alt_rows.append(dict(kind="medium_stem", stem=stem, mapped_exchange=rid, mapped_hg_metabolite=hid, mapped_hg_name=hi["name"],
                                 alt_hg_metabolite=a[0], alt_name=a[1]["name"], alt_formula=a[1]["formula"], alt_charge=a[1]["charge"],
                                 alt_exchange=",".join(hg_ex[a[0]]), alt_kegg=";".join(a[1]["kegg"]), strong_tags=a[2], weak_tags=a[3]))
    else:
        row["class"] = "NOT_IN_SECOND_ROOT"
        if mi:
            row.update(main_name=mi["name"], main_formula=mi["formula"], main_kegg=";".join(mi["kegg"]))
    rows.append(row)

# candidates for the stems missing in the second root. The search keys (KEGG identifiers, name patterns, formulas) are
# search terms chosen for each nutrient; a hit is reported with every identifier the file gives for it.
SEARCH = {
    "adocbl": dict(kegg={"C00194"}, pat=r"cobal|cobam|b12|cbl", forms=set()),
    "4abz": dict(kegg={"C00568"}, pat=r"aminobenzo|\bpaba\b|pabc", forms=set()),
    "no3": dict(kegg={"C00244"}, pat=r"nitrat|\bno3", forms={"NO3"}),
    "zn2": dict(kegg={"C00038"}, pat=r"zinc|\bzn", forms={"Zn"}),
    "mn2": dict(kegg={"C00034"}, pat=r"mangan|\bmn2|\bmn$", forms={"Mn"}),
}
for stem, cfg in SEARCH.items():
    for mid, t in M.hgt.items():
        m = M.sm[mid]
        hit = []
        if set(split_ids(t["metKEGGID"])) & cfg["kegg"]: hit.append("kegg")
        if m.get("formula") in cfg["forms"]: hit.append("formula")
        if re.search(cfg["pat"], m["name"].lower()): hit.append("name")
        if re.search(cfg["pat"], t["metBiGGID"].lower()): hit.append("metBiGGID")
        if not hit: continue
        exs = [r for r in M.ex_by_met.get(mid, [])]
        cand_rows.append(dict(stem=stem, hg_metabolite=mid, compartment=m["compartment"], name=m["name"], formula=m.get("formula", ""),
                              charge=m.get("charge", ""), kegg=t["metKEGGID"], chebi=t["metChEBIID"], metBiGGID=t["metBiGGID"],
                              hit_by=",".join(hit), exchange_in_second_root=";".join(r["pipeline_id"] for r in exs),
                              humangem_reaction=";".join(r["humangem_reaction"] for r in exs), rule_basis=";".join(r["basis"] for r in exs),
                              main_model_has_metabolite=str("%s_e" % stem in M.mm), main_model_has_exchange=str("EX_%s_e" % stem in M.mr)))
    # main-model information on the stem
    for mid in (stem + "_e", stem + "_c"):
        pass

# the six added inputs
for a in added_rows:
    rid = a["exchange_id"]; stem = re.sub(r"^EX_|_e$", "", rid)
    row = dict(kind="added_input", list="added", stem=stem, category="added_by_second_root_build", exchange_id=rid,
               in_main_model=str(rid in M.mr), in_second_root=str(rid in M.sr),
               opened_by_main_list=str(any(s == stem for s, c, l in medium_stems) and rid in M.mr),
               added_label=a["metabolite_name"], added_hg_reaction=a["humangem_reaction"], added_note=a["basis"])
    hid = add_ex_checks(row, rid)
    row["hg_metabolite"] = hid or ""
    hi = M.hg_info(hid) if hid else None
    # reaction in map must equal the one named in medium_additions
    row["ex_reaction_matches_additions"] = str(row.get("ex_reaction_in_map") == a["humangem_reaction"])
    if rid in M.mr:
        mr = M.mr[rid]
        row.update(main_ex_metabolite=";".join(mr["metabolites"]), main_ex_lower=mr["lower_bound"], main_ex_upper=mr["upper_bound"])
    mid, how = main_counterpart(stem, hi)
    row["main_metabolite"] = mid or ""; row["main_counterpart_basis"] = how
    mi = M.main_info(mid) if mid else None
    fill(row, mi, hi)
    sm = search_main(hi) if hi else []
    row["main_exchanges_same_compound_candidates"] = ";".join(
        "%s(%s|%s|%s|%s)" % (",".join(e), i["name"], i["formula"], t, classify(i, hi)[0] if hi else "") for m_, e, i, t in sm)
    al = alternatives(M, hg_ex, hg_e_info, hid, mi, "") if hid else []
    strong = [x for x in al if x[2]]
    row["n_alt_strong"] = len(strong); row["n_alt_weak"] = len(al) - len(strong)
    row["alt_strong"] = ";".join("%s(%s|%s)" % (x[0], ",".join(hg_ex[x[0]]), x[1]["name"]) for x in strong)
    for x in al:
        alt_rows.append(dict(kind="added_input", stem=stem, mapped_exchange=rid, mapped_hg_metabolite=hid, mapped_hg_name=hi["name"],
                             alt_hg_metabolite=x[0], alt_name=x[1]["name"], alt_formula=x[1]["formula"], alt_charge=x[1]["charge"],
                             alt_exchange=",".join(hg_ex[x[0]]), alt_kegg=";".join(x[1]["kegg"]), strong_tags=x[2], weak_tags=x[3]))
    rows.append(row)


def wr(path, rs):
    fields = []
    for r in rs:
        for k in r:
            if k not in fields: fields.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader(); w.writerows(rs)

wr(out / "medium_audit.tsv", rows)
wr(out / "medium_missing_candidates.tsv", cand_rows)
wr(out / "medium_alternatives.tsv", alt_rows)
print(collections.Counter((r["kind"], r["class"]) for r in rows))
print("opened_by_main_list", sum(r["opened_by_main_list"] == "True" for r in rows if r["kind"] == "medium_stem"))
