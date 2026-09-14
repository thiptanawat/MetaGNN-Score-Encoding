"""MCDB 131 culture medium for the Copeland et al. 2023 lung-fibroblast experiment.

This mirrors the structure of the study's own medium module so that the only thing that
changes between the historical cancer-line analysis and this independent evaluation is the
declared culture environment, not the way a culture environment is represented.

Reported medium, Copeland et al. 2023, Materials and methods, "Metabolic flux protocol":
    "MCDB131 medium lacking glucose, glutamine, and phenol red (genDEPOT) which was
     supplemented with 2% dialyzed fetal bovine serum (Mediatech) and naturally labeled
     glucose and glutamine ... For LFs, glucose was supplemented at 8 mM and glutamine was
     supplemented at 1 mM."
    Culture was in 35 mm dishes, 25,000 cells seeded on Day -1, either a standard incubator
    (21 percent oxygen) or a glovebox at 0.5 percent oxygen, both at 5 percent CO2.

Primary sources for the formulation:
    Thermo Fisher / Gibco media formulation 85, MCDB 131 without glutamine:
      https://www.thermofisher.com/us/en/home/technical-resources/media-formulation.85.html
    Caisson Labs MBL02 formulation sheet:
      https://caissonlabs.com/wp-content/uploads/MBL02-formulation.pdf
    HiMedia AL133A technical datasheet:
      https://www.himedialabs.com/media/TD/AL133A.pdf
    Copeland et al. 2023, eLife 12:e82597, https://doi.org/10.7554/eLife.82597

Two decisions that change what this medium means, both declared rather than assumed silently:

  Dialyzed serum.  Dialysis against a 3.5 to 10 kDa membrane removes free small molecules, so
  2 percent dialyzed FBS is represented as contributing no free amino acids, sugars, vitamins
  or ions. What it does retain is macromolecular, so lipids carried on lipoproteins are given a
  small assumed capacity under their own category. This is the main compositional difference
  from the historical RPMI-plus-5-percent-FBS environment, where serum small molecules were
  assumed present.

  Oxygen.  Atmospheric percentage is not a measured uptake bound, so the same assumed oxygen
  capacity is used for both oxygen conditions in the primary analysis and the drug contrasts
  are taken within an oxygen level, where that assumption cancels. A declared sensitivity
  varies the capacity instead of pretending 0.5 percent maps onto a known rate.

Components of MCDB 131 that the pinned Recon3D does not represent are listed in the missing
ledger and are never substituted with a chemically different species: vitamin B12
(cyanocobalamin is absent from the model; aquacobalamin is a different vitamer), calcium and
zinc (present as metabolites but with no exchange reaction), and magnesium, copper, manganese,
nickel, selenite, silicate, vanadate and molybdate (absent from the model entirely).
"""
# (stem, category). Category records availability evidence, never a measured flux bound.
BASAL = [
 # amino acids; MCDB 131 supplies free L-cysteine, not the cystine used in RPMI-1640
 ("gly", "documented_basal"), ("ala__L", "documented_basal"), ("arg__L", "documented_basal"),
 ("asn__L", "documented_basal"), ("asp__L", "documented_basal"), ("cys__L", "documented_basal"),
 ("glu__L", "documented_basal"), ("his__L", "documented_basal"), ("ile__L", "documented_basal"),
 ("leu__L", "documented_basal"), ("lys__L", "documented_basal"), ("met__L", "documented_basal"),
 ("phe__L", "documented_basal"), ("pro__L", "documented_basal"), ("ser__L", "documented_basal"),
 ("thr__L", "documented_basal"), ("trp__L", "documented_basal"), ("tyr__L", "documented_basal"),
 ("val__L", "documented_basal"),
 # vitamins and cofactors
 ("btn", "documented_basal"), ("chol", "documented_basal"), ("pnto__R", "documented_basal"),
 ("5fthf", "documented_basal"), ("ncam", "documented_basal"), ("pydxn", "documented_basal"),
 ("ribflv", "documented_basal"), ("thm", "documented_basal"), ("inost", "documented_basal"),
 ("lipoate", "documented_basal"),
 # inorganic ions with an exchange in the pinned model
 ("fe2", "documented_basal"), ("k", "documented_basal"), ("na1", "documented_basal"),
 ("pi", "documented_basal"), ("so4", "documented_basal"), ("cl", "documented_basal"),
 ("hco3", "documented_basal"),
 # other defined organics
 ("ade", "documented_basal"), ("ptrc", "documented_basal"), ("pyr", "documented_basal"),
 ("thymd", "documented_basal"),
 # formulation components with no chemically exact exchange in the pinned model
 ("__b12__", "unrepresented_source_b12"),
 ("ca2", "unrepresented_trace"), ("zn2", "unrepresented_trace"), ("mg2", "unrepresented_trace"),
 ("cu2", "unrepresented_trace"), ("mn2", "unrepresented_trace"), ("ni2", "unrepresented_trace"),
 ("sel", "unrepresented_trace"), ("sio3", "unrepresented_trace"), ("vanad", "unrepresented_trace"),
 ("moo4", "unrepresented_trace"),
]
CULTURE = [
 ("glc__D", "documented_culture"),          # supplemented at 8 mM for lung fibroblasts
 ("gln__L", "documented_culture"),          # supplemented at 1 mM for lung fibroblasts
 ("o2", "inferred_aerobic_availability"),   # assumed capacity, identical across oxygen conditions
 ("h2o", "assumed_solvent_exchange"),
 ("h", "assumed_proton_exchange"),
]
SERUM = [  # 2 percent DIALYZED FBS: free small molecules removed, lipoprotein-bound lipids retained
 ("chsterol", "assumed_lipoprotein_bound"), ("ocdcea", "assumed_lipoprotein_bound"),
 ("hdca", "assumed_lipoprotein_bound"), ("ocdca", "assumed_lipoprotein_bound"),
 ("lnlc", "assumed_lipoprotein_bound"), ("arachd", "assumed_lipoprotein_bound"),
]
UPTAKE = {"documented_basal": 10.0, "documented_culture": 10.0,
          "assumed_lipoprotein_bound": 1.0, "inferred_aerobic_availability": 10.0,
          "assumed_solvent_exchange": 10.0, "assumed_proton_exchange": 10.0}

AVAILABILITY_BASIS = {
 "documented_basal": "supplier MCDB 131 formulation tables (Gibco 85, Caisson MBL02, HiMedia AL133A); the exact genDEPOT lot is not published",
 "documented_culture": "Copeland et al. 2023 methods state 8 mM glucose and 1 mM L-glutamine for lung fibroblasts",
 "assumed_lipoprotein_bound": "2 percent dialyzed FBS retains lipoproteins; the lipid identities and capacities assigned to them are assumptions",
 "inferred_aerobic_availability": "oxygen capacity assumed and held identical across oxygen conditions; atmospheric percentage is not a measured uptake bound",
 "assumed_solvent_exchange": "aqueous-solvent modeling allowance, not measured water transport",
 "assumed_proton_exchange": "proton-balance modeling allowance, not a pH constraint",
 "unrepresented_source_b12": "cyanocobalamin is absent from the pinned model; aquacobalamin is a different B12 vitamer and is not substituted",
 "unrepresented_trace": "formulation component with no chemically exact exchange reaction in the pinned model",
}


def apply_medium(model, uptake=None):
    """Apply declared availability assumptions; return opened, missing and a bound ledger."""
    up = dict(UPTAKE); up.update(uptake or {})
    ledger = []
    for r in model.reactions:
        if r.id.startswith(("DM_", "SK_")):
            if r.lower_bound < 0.0:        # an intracellular free source would bypass the medium
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
        if cat in ("unrepresented_source_b12", "unrepresented_trace"):
            if rid in model.reactions:
                raise ValueError("The pinned model changed: %s now exists and its chemical "
                                 "identity must be resolved before it is left closed" % rid)
            missing.append({"exchange": rid, "category": cat,
                            "availability_basis": AVAILABILITY_BASIS[cat]})
            continue
        if rid in model.reactions:
            lim = up[cat]
            model.reactions.get_by_id(rid).lower_bound = -lim
            opened[rid] = {"lower_bound": -lim, "category": cat,
                           "availability_basis": AVAILABILITY_BASIS[cat],
                           "bound_basis": "assumed transport capacity, not converted from the medium concentration"}
            ledger.append((rid, "lower_bound", 0.0, -lim, cat))
        else:
            missing.append({"exchange": rid, "category": cat,
                            "note": "no such exchange in the pinned model",
                            "availability_basis": AVAILABILITY_BASIS[cat]})
    return opened, missing, ledger
