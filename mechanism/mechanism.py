"""Mechanism analyses for Paper 3, 14 September 2026.

Two questions the constraint-stage pilot left open, both answerable from saved outputs
with no new optimization:

  Q1  Within-profile ranges contract sharply once the primary cost cap is imposed. Does the
      contracted interval sit at the SAME location in every profile? If it does, the objective
      collapses every profile onto one point of the shared admissible set and across-profile
      constancy follows. If the locations differ, constancy has some other cause.

  Q2  Do the encodings actually select different flux vectors, and if so, does that difference
      survive the projection onto the reported exchange coordinates?

Everything here reads frozen artefacts. It reads no CORE outcome value and solves no LP.

Conventions fixed before any output was inspected:
  * A coordinate is numerically FIXED when its range width is <= WIDTH_TOL.
  * Interval LOCATION is the midpoint (lo + hi) / 2.
  * Across-profile spread of a quantity is max - min over profiles.
  * Exploration of a target = across-profile location spread / shared B0 width, defined only
    when the B0 width exceeds WIDTH_TOL.
  * Regret of vector v under cost c with optimum E* is (c.v - E*) / E*, defined for E* > 0.
"""
import numpy as np

WIDTH_TOL = 1e-5          # the study's own numerical-fixedness threshold
REGRET_TOL = 1e-9         # relative slack treated as exact optimality


def width(lo, hi):
    return float(hi) - float(lo)


def location(lo, hi):
    return 0.5 * (float(lo) + float(hi))


def is_fixed(lo, hi, tol=WIDTH_TOL):
    return width(lo, hi) <= tol


def spread(values):
    """Across-profile spread. Empty input has no spread."""
    v = [float(x) for x in values]
    if not v:
        return float("nan")
    return max(v) - min(v)


def exploration(locations, b0_lo, b0_hi, tol=WIDTH_TOL):
    """Fraction of the shared admissible width that the profiles actually move across.

    Returns nan when the shared width is itself numerically fixed, because the ratio is then
    undefined rather than zero.
    """
    w = width(b0_lo, b0_hi)
    if w <= tol:
        return float("nan")
    return spread(locations) / w


def linear_cost(flux, costs):
    """The study's primary objective: sum of cost-weighted absolute flux."""
    return float(np.dot(np.asarray(costs, dtype=float), np.abs(np.asarray(flux, dtype=float))))


def regret(flux, costs, optimum, tol=REGRET_TOL):
    """Relative excess of this vector over the optimum of the cost it is scored against."""
    if not np.isfinite(optimum) or optimum <= 0:
        return float("nan")
    r = (linear_cost(flux, costs) - float(optimum)) / float(optimum)
    return 0.0 if abs(r) <= tol else r


def relative_l1(a, b):
    """Relative L1 distance, symmetric, zero for identical vectors, nan when both are zero."""
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    denom = np.abs(a).sum() + np.abs(b).sum()
    if denom == 0.0:
        return float("nan")
    return float(2.0 * np.abs(a - b).sum() / denom)
