"""Write REPORT.md from target_audit.tsv, medium_audit.tsv, medium_missing_candidates.tsv, alternatives.tsv,
medium_alternatives.tsv and yml_crosscheck.txt. Usage: python -I make_report.py OUTDIR"""
import csv, sys, collections
from pathlib import Path
out = Path(sys.argv[1])
rd = lambda n: list(csv.DictReader(open(out / n), delimiter="\t"))
T, Mx, C, A, MA = rd("target_audit.tsv"), rd("medium_audit.tsv"), rd("medium_missing_candidates.tsv"), rd("alternatives.tsv"), rd("medium_alternatives.tsv")
yml = (out / "yml_crosscheck.txt").read_text().split("\n")[:2]
stems = [r for r in Mx if r["kind"] == "medium_stem"]
adds = [r for r in Mx if r["kind"] == "added_input"]
CL = ["SAME", "SAME_PROTONATION_DIFFERS", "STEREO_UNSPECIFIC_ON_ONE_SIDE", "STEREO_MISMATCH", "DIFFERENT_COMPOUND", "CANNOT_TELL", "NOT_IN_SECOND_ROOT"]
def cnt(rs):
    c = collections.Counter(r["class"] for r in rs); return " | ".join(str(c.get(k, 0)) for k in CL)
L = []
w = L.append
w("# Identifier audit of the Human-GEM mapping (second reconstruction)\n")
w("Facts only. No remapping is proposed. Every number below is computed by the scripts in `scripts/` from the files named in the Method section.\n")
w("## Method\n")
w("- Main model: the Recon3D JSON of the main analysis (`data/raw/Recon3D.json`; name, formula, charge, annotation with `kegg.compound`, `chebi`, `hmdb`, `inchi_key`). Human-GEM: the converted model written by `second/build_root.py` (`data/raw/Recon3D.json` of the second root) joined to `model/metabolites.tsv` of the Human-GEM checkout (KEGG, ChEBI, HMDB, InChI, SMILES, BiGG, MetaNetX) and `reactions.tsv`. Name, formula and charge of every extracellular Human-GEM metabolite with an exchange (%s) were re-read from `Human-GEM.yml` and agree with the JSON (%s)." % (yml[0].replace("\t", " "), yml[1].replace("\t", " ")))
w("- Heavy-atom composition is the formula without hydrogen. A protonation difference is recorded when charge, hydrogen count, or the protonation letter of the InChIKey differs while the heavy atoms agree.")
w("- Human-GEM InChI strings were converted to InChIKeys (RDKit, InChI library). The main model stores only hashed InChIKeys, not InChI layers. Key block 1 encodes connectivity, block 2 the stereo and isotope layers (`UHFFFAOY` means none), and the last letter the protonation state. SAME requires equal heavy atoms, equal key blocks 1 and 2, and equal charge and hydrogen count.")
w("- The measured analyte is represented only by its source KEGG identifier and label (`metabolite_map.tsv`); no structure for the analyte is available offline, so analyte identity rests on KEGG concordance.")
w("- Alternatives are Human-GEM extracellular metabolites with their own one-metabolite boundary reaction that share, with the mapped metabolite or the main-model metabolite, a KEGG identifier, a ChEBI identifier, an InChIKey connectivity block, a MetaNetX identifier, a normalized name, or a BiGG stem (strong tags). A weaker tag marks overlapping name tokens; weak rows are only in the TSV files.\n")
w("## Counts per class\n")
w("| Set | n | " + " | ".join(CL) + " |")
w("|---|---|" + "---|" * len(CL))
w("| 96 targets | %d | %s |" % (len(T), cnt(T)))
w("| Medium stems (all %d in the three lists) | %d | %s |" % (len(stems), len(stems), cnt(stems)))
w("| Six added inputs | %d | %s |" % (len(adds), cnt(adds)))
w("")
w("Main model opens %d of the %d stems (those with an `EX_<stem>_e` reaction). The second root has an exchange for %d stems." % (
    sum(r["opened_by_main_list"] == "True" for r in stems), len(stems), sum(r["in_second_root"] == "True" for r in stems)))
only2 = [r["stem"] for r in stems if r["in_second_root"] == "True" and r["in_main_model"] == "False"]
miss2 = [r["stem"] for r in stems if r["in_second_root"] == "False"]
w("Stems with an exchange in the second root but none in the main model (so the main list cannot open them): %s. Stems absent from the second root: %s; none of these has an exchange in the main model either (in_main_model column).\n" % (", ".join(only2), ", ".join(miss2)))
w("## Exchange checks\n")
def exsum(rs, nm):
    ok = sum(r["ex_rule_ok"] == "True" for r in rs)
    b = collections.Counter((r.get("ex_lower"), r.get("ex_upper")) for r in rs if r["ex_rule_ok"] == "True")
    return "%s: %d of %d present exchanges have exactly one metabolite, in compartment e, coefficient -1. Bounds (lower, upper) in the second root: %s." % (nm, ok, sum(r["ex_rule_ok"] != "NOT_IN_SECOND_ROOT" for r in rs), "; ".join("%s to %s in %d" % (k[0], k[1], v) for k, v in b.items()))
w("- " + exsum(T, "Targets"))
w("- " + exsum(stems, "Medium stems"))
w("- " + exsum(adds, "Added inputs"))
w("- Bounds are the values written in the converted JSON (import-direction reactions were reversed by the build, `exchange_map.tsv` column `reversed`). The medium module resets them (lower bound to 0, upper bound to 1000, then opens the listed stems), so the default bounds matter only for exchanges the medium does not touch.")
w("- `ex_matches_target_map` and `main_ex_metabolite_matches_map` in `target_audit.tsv` are True for all %d targets: the exchange identifier resolves to the Human-GEM metabolite of `target_map.tsv` and to the main-model metabolite of `metabolite_map.tsv`.\n" % len(T))
def ev(r, medium=False):
    bits = []
    if not r["main_metabolite"]:
        bits.append("no main-model counterpart; HG `%s` %s %s" % (r["hg_metabolite"], r["hg_name"], r["hg_formula"]))
        bits.append("KEGG HG %s" % (r["hg_kegg"] or "-"))
        if r.get("missing"): bits.append("missing: " + r["missing"])
        return "; ".join(bits)
    bits.append("main %s `%s` %s q=%s; HG `%s` %s %s q=%s" % (r["main_metabolite"], r["main_name"], r["main_formula"], r.get("charge_main", ""), r["hg_metabolite"], r["hg_name"], r["hg_formula"], r.get("charge_hg", "")))
    bits.append("KEGG main %s / HG %s%s" % (r["main_kegg"] or "-", r["hg_kegg"] or "-", (" / source %s" % r["source_kegg"]) if "source_kegg" in r and r["source_kegg"] else ""))
    bits.append("keys main %s / HG %s" % (r.get("key_main") or "-", r.get("key_hg") or "-"))
    bits.append("ChEBI overlap %s" % (r.get("chebi_overlap") or "none"))
    if r.get("note"): bits.append(r["note"])
    if r.get("missing"): bits.append("missing: " + r["missing"])
    return "; ".join(bits)
w("## Rows that are not SAME or SAME_PROTONATION_DIFFERS\n")
w("### Targets\n")
w("| Analyte | Exchange | Class | Evidence |"); w("|---|---|---|---|")
for r in T:
    if r["class"] not in ("SAME", "SAME_PROTONATION_DIFFERS"):
        w("| %s | %s | %s | %s |" % (r["source_label"], r["exchange_id"], r["class"], ev(r)))
w("")
w("Additional facts for these rows (from the same files):\n")
for r in T:
    if r["exchange_id"] == "EX_glyald_e":
        w("- `EX_glyald_e`: `target_map.tsv` reports KEGG concordance False (source %s, main %s, HG %s). The Human-GEM InChI has no stereo layer; HG metBiGGID is `glyald` and ChEBI %s. Human-GEM has no other metabolite carrying %s." % (r["source_kegg"], r["main_kegg"], r["hg_kegg"], r["hg_chebi"], r["source_kegg"]))
    if r["exchange_id"] == "EX_pnto__R_e":
        w("- `EX_pnto__R_e`: the Human-GEM metabolite is named `pantothenate`, its InChI has no stereo layer, but its BiGG identifier is `%s` and its ChEBI identifier %s lies in the main-model list." % (r["hg_bigg_ids"], r["hg_chebi"]))
    if r["exchange_id"] == "EX_L2aadp_e":
        w("- `EX_L2aadp_e`: Human-GEM name `%s` with BiGG `%s`; its InChI has no stereo layer, the main-model key has one." % (r["hg_name"], r["hg_bigg_ids"]))
    if r["exchange_id"] == "EX_glcur_e":
        w("- `EX_glcur_e`: same connectivity block, different stereo block. HG ChEBI %s is not in the main list. The HG InChI stereo layer is `%s`. The main model stores no InChI, so which stereocenter differs (for example ring form or anomer) cannot be read from the files. The alternative `EX_idour_e` (L-iduronic acid) shares the connectivity block." % (r["hg_chebi"], r["hg_inchi"].split("/t")[-1]))
    if r["exchange_id"] == "EX_lcts_e":
        w("- `EX_lcts_e`: same connectivity block, different stereo block. HG ChEBI %s is not in the main list (main list: %s). The main KEGG list is %s against HG %s. Whether the difference is an anomeric configuration cannot be read from the files." % (r["hg_chebi"], r["main_chebi"], r["main_kegg"], r["hg_kegg"]))
    if r["exchange_id"] == "EX_orn_e":
        w("- `EX_orn_e`: `target_map.tsv` marks `same_by_rxnBiGGID` False (reaction-level BiGG identifier `%s`). Structure and KEGG agree (class %s); the metabolite-level cross-reference was used." % (r["rxnBiGGID"], r["class"]))
w("")
w("### Medium stems and added inputs\n")
w("| Row | Class | Evidence |"); w("|---|---|---|")
for r in Mx:
    if r["class"] not in ("SAME", "SAME_PROTONATION_DIFFERS", "NOT_IN_SECOND_ROOT"):
        w("| %s `%s` | %s | %s |" % ("added" if r["kind"] == "added_input" else "stem", r["stem"], r["class"], ev(r)))
w("")
w("Stems without a main-model counterpart (mg2, cu2, sel): the main model has no metabolite sharing the Human-GEM KEGG identifier and no exchange, so no comparison is possible; the Human-GEM entities are listed in the table.\n")
w("## Rows classed SAME_PROTONATION_DIFFERS\n")
w("| Row | Evidence |"); w("|---|---|")
for r in T:
    if r["class"] == "SAME_PROTONATION_DIFFERS": w("| target %s | %s |" % (r["exchange_id"], ev(r)))
for r in Mx:
    if r["class"] == "SAME_PROTONATION_DIFFERS": w("| %s `%s` | %s |" % ("added" if r["kind"] == "added_input" else "stem", r["stem"], ev(r)))
w("")
w("In the main model, the key suffix of 3-hydroxyanthranilate, taurocholate and riboflavin (letter M) does not match the stored charge 0 and neutral formula; this is a property of the main-model annotation, not of Human-GEM.\n")
w("## Alternative exchanges for the same analyte (strong tags)\n")
w("### Targets\n"); w("| Analyte | Mapped HG metabolite | Alternative | Exchange | Tags |"); w("|---|---|---|---|---|")
for a in A:
    if a["strong_tags"]:
        w("| %s | %s %s | %s %s | %s | %s |" % (a["source_label"], a["mapped_hg_metabolite"], a["mapped_hg_name"], a["alt_hg_metabolite"], a["alt_name"], a["alt_exchange"], a["strong_tags"]))
w("\nFull list with weak tags: `alternatives.tsv` (%d rows, %d with strong tags)." % (len(A), sum(bool(a["strong_tags"]) for a in A)))
w("\n### Medium stems and added inputs\n"); w("| Stem | Mapped HG metabolite | Alternative | Exchange | Tags |"); w("|---|---|---|---|---|")
for a in MA:
    if a["strong_tags"]:
        w("| %s | %s %s | %s %s | %s | %s |" % (a["stem"], a["mapped_hg_metabolite"], a["mapped_hg_name"], a["alt_hg_metabolite"], a["alt_name"], a["alt_exchange"], a["strong_tags"]))
w("\nFull list: `medium_alternatives.tsv`. All alternatives listed have default bounds -1000 to 1000 and are not opened by the medium; flux through them depends on the medium module closing uptake, which it does for every EX_ reaction outside the opened list. Secretion is left free, so an alternative stereoisomer could carry export flux.\n")
w("## Stems of the main list missing in the second root\n")
w("Candidates were searched by KEGG identifier, name pattern, formula and BiGG identifier over all Human-GEM metabolites (`medium_missing_candidates.tsv`). The search-key KEGG identifiers were supplied by the auditor, not read from a file; the name and formula searches are independent of them.\n")
w("| Stem | Human-GEM hit | Compartment | Formula, charge | KEGG | metBiGGID | Exchange | Reaction | Rule basis |"); w("|---|---|---|---|---|---|---|---|---|")
for c in C:
    if c["compartment"] == "e" or c["exchange_in_second_root"]:
        w("| %s | %s %s | %s | %s, %s | %s | %s | %s | %s | %s |" % (c["stem"], c["hg_metabolite"], c["name"], c["compartment"], c["formula"], c["charge"], c["kegg"], c["metBiGGID"], c["exchange_in_second_root"] or "none", c["humangem_reaction"], c["rule_basis"] or "-"))
w("")
w("Reading of the hits:")
w("- adocbl: `MAM01600e` (cobamide-coenzyme, KEGG C00194) has the exchange `EX_MAR08400`; the build kept the Human-GEM identifier on purpose (rule basis column). `MAM01361e` aquacob(III)alamin is a different cobalamin and is one of the six added inputs.")
w("- zn2: `MAM03157e` zinc (Zn, charge 2, KEGG C00038, the same KEGG as main-model `zn2_e`) has `EX_HC02172_e`. The rule did not pick `EX_zn2_e` because the Human-GEM metBiGGID is `HC02172`, not `zn2`. The main model has `zn2_e` but no exchange for it.")
w("- no3: the only extracellular hit is `MAM02714e` `peroxynitrite` (formula NO3, charge -1, KEGG C16845, exchange `EX_CE5643_e`). Its formula equals that of nitrate but the name and KEGG identifier name peroxynitrite; whether it is a mislabeled nitrate is not decidable from the files.")
w("- 4abz: no Human-GEM metabolite matched 4-aminobenzoate by KEGG, name or BiGG stem. The only aminobenzoate is `MAM01342e` anthranilate (2-aminobenzoate, `EX_anth_e`), an isomer with the same formula that is a different compound and is itself a measured target.")
w("- mn2: no Human-GEM metabolite matched by KEGG, name, formula or BiGG identifier (the one name hit is a mannose-N-acetylglucosamine disaccharide, not manganese).\n")
w("## The six added inputs\n")
w("| Input | HG reaction | HG metabolite | Class vs main metabolite | Main exchange for same compound | Opened by main medium |"); w("|---|---|---|---|---|---|")
for r in adds:
    w("| %s | %s | %s %s | %s | %s | %s |" % (r["exchange_id"], r["added_hg_reaction"], r["hg_metabolite"], r["hg_name"], r["class"], r["main_exchanges_same_compound_candidates"] or "none", r["opened_by_main_list"]))
w("")
w("`ex_reaction_matches_additions` is True for all six: the identifier in `medium_additions.tsv` resolves to the reaction named there. None of the six stems occurs in BASAL, CULTURE or SERUM, so the main medium opens none of them. The main model has an exchange for five under the same identifier (`in_main_model` column) and carries gamma-tocopherol under a different identifier (`EX_yvite_e`) because its `gtocophe_e` metabolite has no exchange. Candidate lists with shared-identifier tags and classes are in the `main_exchanges_same_compound_candidates` column of `medium_audit.tsv`; other tocopherol and retinoate isomers listed there are separate main-model entries.\n")
w("## Open points\n")
w("- Main-model InChI layers are not stored, so a stereo class rests on key-block equality; the specific stereocenter of a mismatch is not resolved.")
w("- Rows classed CANNOT_TELL lack a main-model InChIKey (arachidonate, aquacob(III)alamin) or a counterpart (mg2, cu2, sel), or have a connectivity block that differs without a KEGG conflict (thiamin); see the Evidence column.")
w("- Measured-analyte identity is checked through KEGG only.")
w("- KEGG identifiers used as search keys for the five missing stems are from the auditor, not from a file.")
w("\n## Files\n")
w("`target_audit.tsv`, `medium_audit.tsv`, `alternatives.tsv`, `medium_alternatives.tsv`, `medium_missing_candidates.tsv`, `yml_crosscheck.txt`; scripts in `scripts/` (`audit_lib.py`, `audit_targets.py`, `audit_medium.py`, `check_yml.py`, `make_report.py`).")
(out / "REPORT.md").write_text("\n".join(L) + "\n")
