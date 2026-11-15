# Plan: the quantities of the second reconstruction without the profiles that the second solver completed

Fixed on 15 November 2026, before the computation it describes. It adds to `second/SECOND_RECONSTRUCTION_PLAN_2026-11-14.md` and changes nothing in it.

## 1. Question

In the reserved run on Human-GEM the solver of the main analysis (GLPK) stopped on an infeasible problem in 4 of the 283 arms, and the second solver (HiGHS) completed them, as Section 3 of that plan provides. The points of those 4 arms entered every comparison beside points of the first solver. The check of the two solvers (step F of that plan) compares them on arms that both solved, so it does not show what the 4 arms contribute. Do the reported quantities depend on them?

## 2. What is removed

The 4 arms are OVCAR-4 under the identity encoding, HCT-15 and IGROV1 under the power encoding with exponent 1.25, and MDA-MB-435 under the midrank encoding (`results/second/reserved_v3/second_solver_arms.json`). The evaluation of the main analysis gives no contrast when an arm of a profile is missing. The 4 profiles are therefore removed whole, under every encoding, which leaves 43 of the 47 reserved profiles. Nothing else changes: the targets, the tolerances, the pair rule, the seed and the number of bootstrap resamples are those of the locked analysis plan of the second root.

## 3. What is computed

`second/without_second_solver.py`, committed with this plan, computes the following for the 43 profiles from the saved records of the reserved run. No linear program is solved, and no file of the reserved run is changed.

| Statement | Quantity | How |
|---|---|---|
| Encoding | Concordance of every encoding on the 52 primary targets; the contrasts midrank minus the mean of the power encodings and midrank minus identity, with their 95% bootstrap intervals over origins | `src/evaluation.py` (`evaluate_panel`), unchanged |
| Encoding | Predictions changed by the encoding and between-origin orderings reversed | `src/transform_sensitivity.py` (`summarize`), unchanged |
| Attribution | Coordinates that the network leaves open and coordinates first fixed by the primary cost cap, under every encoding | the rows of `results/second/mechanism/reserved_range/summary.json` of the 43 profiles |
| Attribution | Targets constant across the profiles, and how many of them take the value of the uniform solution | the rule of `second/summarize.py` |
| Allowance | Mean number of fixed coordinates and share of between-origin orderings certified by the interval rule, under the near-exact cap and under the caps of 1%, 5% and 20%, for the identity and the midrank encoding | the job records of `results/second/closure/cap_sweep/` and the rule of `closure/cap_concordance.py` |

Before anything is removed, the script computes every quantity for all 47 profiles with the same code and compares it with the deposited value. It stops if one differs.

## 4. Reading

The three readings of Section 6 of the plan of the second reconstruction are applied to the 43 profiles as they were applied to the 47. The manuscript reports for each statement whether it holds without the 4 profiles, and gives the two concordance contrasts with their intervals beside those of all profiles, whatever they are. A bootstrap over fewer origins has wider intervals; the reading of the encoding statement is the one of that plan and is not adjusted for this.

## 5. What was known when this plan was fixed

Every result of the second reconstruction for all 47 profiles, which release v1.3.0 holds. The script had been run once with `--check-only`, which computes the quantities for all 47 profiles, compares them with the deposited values (all equal) and removes nothing. No quantity without the 4 profiles had been computed.

## 6. What is not recomputed

The agreement of the two solvers (step F), which removes no profile. The range calculations and the cap sweep themselves: their saved records are read as they are, and those of the 4 profiles are left out. The uniform-cost control, which does not depend on a profile. The interval-position analysis (`mechanism/range_location.py`), whose statistics over the profiles are not computed again.

## 7. Deviations

None at the time the plan was fixed.
