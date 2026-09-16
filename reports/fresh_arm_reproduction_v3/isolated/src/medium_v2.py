"""Source-informed RPMI/serum approximation for the pinned generic Recon3D.

Availability provenance and numeric uptake assumptions are separate. CORE's
supplement p. 2 reports GIBCO RPMI-1640, 2 mM glutamine, 5% FBS and 5% CO2/95% air,
but no catalog number or lot. The current GIBCO 11875 formulation is a representative
reference, not proof of the exact 2012 preparation. No measured CORE exchange rate,
sign or cell-specific quantity is used. The medium is identical across contexts.

Primary sources:
  CORE culture: https://pmc.ncbi.nlm.nih.gov/articles/PMC3526189/ (supplement p. 2)
  Representative GIBCO formulation:
    https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.114.html
  Bicarbonate buffering, product information p. 1:
    https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0018935_RPMI_1640Medium_IFU.pdf
  Form-explicit representative RPMI vitamin B12 (a different manufacturer):
    https://bioscience.lonza.com/lonza_bs/US/en/document/30050
  Adenosylcobalamin identity: https://www.kegg.jp/entry/C00194

The GIBCO reference says vitamin B12 (MW 1355), consistent with cyanocobalamin;
adenosylcobalamin is a different chemical form (MW 1579.58). The inherited adocbl
placeholder has no exchange in the pinned model and supplies nothing. It is retained
only in the missing-ingredient ledger, never described as an exact ingredient match.
Oxygen availability is inferred from air; its uptake capacity is assumed. Water and
proton exchange are solvent/buffering conventions, not measured transport or pH.

Corrections in this revision:
  * L-cystine (EX_cysi__L_e), not L-cysteine. RPMI-1640 lists L-cystine dihydrochloride.
  * trans-4-hydroxy-L-proline added; it is in the formulation and was previously omitted.
  * Boundary classes now have separate, documented policies rather than one blanket rule.

Boundary policy, stated explicitly because it changes what the medium means:
  EX_  extracellular exchange. Uptake closed, then reopened for declared source-informed
       ingredients or explicit environment/serum assumptions.
  DM_  intracellular demand. Left at the pinned model's bounds; these drain a metabolite and do not
       supply one. Three of the 145 carry a negative lower bound and are closed on the supply side
       for the same reason as the sinks.
  SK_  intracellular sink. 92 of 101 carry a lower bound of -1000 in the pinned model, so they can
       SUPPLY metabolite from nothing. Left open they bypass the medium entirely: biomass then
       reaches 371 rather than 3.65. The supply direction is therefore closed and the drain
       direction left free. This is a modeling decision, recorded here and in the change ledger,
       not a silent repair.
"""
# (stem, category) — category records availability evidence, never a measured flux bound.
BASAL = [  # components supported by the representative GIBCO formulation
 ("arg__L","documented_basal"),("asn__L","documented_basal"),("asp__L","documented_basal"),
 ("cysi__L","documented_basal"),("glu__L","documented_basal"),("gly","documented_basal"),
 ("his__L","documented_basal"),("4hpro_LT","documented_basal"),("ile__L","documented_basal"),
 ("leu__L","documented_basal"),("lys__L","documented_basal"),("met__L","documented_basal"),
 ("phe__L","documented_basal"),("pro__L","documented_basal"),("ser__L","documented_basal"),
 ("thr__L","documented_basal"),("trp__L","documented_basal"),("tyr__L","documented_basal"),
 ("val__L","documented_basal"),
 ("btn","documented_basal"),("pnto__R","documented_basal"),("chol","documented_basal"),
 ("fol","documented_basal"),("inost","documented_basal"),("ncam","documented_basal"),
 ("pydxn","documented_basal"),("ribflv","documented_basal"),("thm","documented_basal"),
 ("adocbl","unrepresented_source_b12"),("4abz","documented_basal"),
 ("glc__D","documented_basal"),("gthrd","documented_basal"),
 ("ca2","documented_basal"),("cl","documented_basal"),("k","documented_basal"),
 ("na1","documented_basal"),("mg2","documented_basal"),("pi","documented_basal"),
 ("so4","documented_basal"),("no3","documented_basal"),
]
CULTURE = [  # distinct direct, inferred, and modeling provenance
 ("gln__L","documented_culture"),           # RPMI-1640 with added L-glutamine
 ("o2","inferred_aerobic_availability"),    # air atmosphere, not a measured uptake rate
 ("h2o","assumed_solvent_exchange"),
 ("h","assumed_proton_exchange"),          # does not impose physiological pH
 ("hco3","documented_basal"),              # representative GIBCO bicarbonate buffer
]
SERUM = [  # 5% fetal bovine serum: present by documentation, composition assumed
 ("chsterol","assumed_serum"),("ocdcea","assumed_serum"),("hdca","assumed_serum"),
 ("ocdca","assumed_serum"),("lnlc","assumed_serum"),("arachd","assumed_serum"),
 ("crn","assumed_serum"),("pyr","assumed_serum"),("lac__L","assumed_serum"),
 ("ac","assumed_serum"),("urea","assumed_serum"),("gthox","assumed_serum"),
 ("fe2","assumed_serum"),("fe3","assumed_serum"),("zn2","assumed_serum"),
 ("cu2","assumed_serum"),("mn2","assumed_serum"),("sel","assumed_serum"),
]
# Uptake capacities are ASSUMED. Medium concentrations constrain availability, not mmol/gDW/h rates.
UPTAKE = {"documented_basal":10.0, "documented_culture":10.0, "assumed_serum":1.0,
          "unrepresented_source_b12":10.0, "inferred_aerobic_availability":10.0,
          "assumed_solvent_exchange":10.0, "assumed_proton_exchange":10.0}

# The legacy documented_basal key is retained for compatibility; its evidence is
# the representative formulation above, not direct ingredient-level CORE records.
AVAILABILITY_BASIS = {
 "documented_basal": "representative current GIBCO RPMI formulation; exact 2012 catalog/lot unknown",
 "documented_culture": "CORE supplement p. 2 directly reports 2 mM L-glutamine",
 "assumed_serum": "CORE reports 5% FBS; the identities and capacities assigned to its constituents are assumptions",
 "unrepresented_source_b12": "source vitamin B12 form not specified; adocbl is not an exact demonstrated match and no exchange is implemented",
 "inferred_aerobic_availability": "oxygen availability inferred from reported 95% air; gas transfer and capacity unmeasured",
 "assumed_solvent_exchange": "aqueous-solvent modeling allowance; not measured water transport",
 "assumed_proton_exchange": "proton-balance modeling allowance; not a pH constraint or measured proton transport",
}

def apply_medium(model, uptake=None):
    """Apply declared availability assumptions; return opened, missing, and bound ledger."""
    up = dict(UPTAKE); up.update(uptake or {})
    ledger = []
    for r in model.reactions:
        if r.id.startswith(("DM_", "SK_")):
            if r.lower_bound < 0.0:                    # an intracellular free source bypasses the medium
                ledger.append((r.id, "lower_bound", r.lower_bound, 0.0,
                               "close the supply direction of an intracellular boundary"))
                r.lower_bound = 0.0
            continue
        if not r.id.startswith("EX_"):
            continue
        if r.lower_bound != 0.0:
            ledger.append((r.id, "lower_bound", r.lower_bound, 0.0, "close uptake by default"))
            r.lower_bound = 0.0
        if r.upper_bound < 1000.0:
            ledger.append((r.id, "upper_bound", r.upper_bound, 1000.0, "secretion left free"))
            r.upper_bound = 1000.0
    opened, missing = {}, []
    for stem, cat in BASAL + CULTURE + SERUM:
        rid = "EX_%s_e" % stem
        if cat == "unrepresented_source_b12":
            if rid in model.reactions:
                raise ValueError("The pinned model changed: resolve vitamin-B12 chemical identity before opening an exchange")
            missing.append({"exchange": rid, "category": cat,
                            "note": "unrepresented source vitamin B12; no adenosylcobalamin substitution is made",
                            "availability_basis": AVAILABILITY_BASIS[cat]})
            continue
        if rid in model.reactions:
            lim = up[cat]
            model.reactions.get_by_id(rid).lower_bound = -lim
            opened[rid] = {"lower_bound": -lim, "category": cat,
                           "availability_basis": AVAILABILITY_BASIS[cat],
                           "bound_basis": "assumed transport capacity, not converted from concentration or CORE"}
            ledger.append((rid, "lower_bound", 0.0, -lim, cat))
        else:
            missing.append({"exchange": rid, "category": cat, "note": "no such exchange in the pinned model",
                            "availability_basis": AVAILABILITY_BASIS[cat]})
    return opened, missing, ledger
