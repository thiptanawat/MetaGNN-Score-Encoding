"""The reporting rule and the measurement-eligibility rule, kept separate from the data.

Sign conventions, stated once and checked by the tests:
  measured   fmol per cell per hour, positive is secretion and negative is uptake
             (compendium R/data.R: "positive fluxes indicate secretion, negative fluxes
             indicate uptake")
  predicted  canonical model units, exchange reactions written as metabolite ->, so a positive
             flux is secretion and a negative flux is uptake
Both conventions agree, so a direction can be compared without any sign flip. Only the
ORDERING is compared; the two unit systems are never equated.

A direction of +1 always means "the first condition is greater than the second" in both.
"""
EPS = 1e-5


def measured_direction(values_a, values_b):
    """+1 when every replicate of A exceeds every replicate of B, -1 for the reverse, else 0.

    Zero means the measurement does not resolve a direction, not that the two are equal.
    """
    if not values_a or not values_b:
        return 0
    if min(values_a) > max(values_b):
        return 1
    if min(values_b) > max(values_a):
        return -1
    return 0


def interval_direction(intervals_a, intervals_b, eps=EPS):
    """Report a direction only when every A interval clears every B interval by eps.

    intervals_* are lists of (lo, hi) over encodings and replicate profiles taken together.
    """
    if not intervals_a or not intervals_b:
        return 0
    if all(lo_a > hi_b + eps for lo_a, _ in intervals_a for _, hi_b in intervals_b):
        return 1
    if all(lo_b > hi_a + eps for lo_b, _ in intervals_b for _, hi_a in intervals_a):
        return -1
    return 0


def shared_encoding_direction(intervals_by_arm_a, intervals_by_arm_b, eps=EPS):
    """The planned rule: every encoding must support the same direction on its own intervals.

    intervals_by_arm_* map an encoding name to the list of (lo, hi) over the replicate profiles
    of one condition. Within each encoding every cross-replicate pair must clear eps; a
    direction is reported only when every encoding in the set reports that same direction.
    Intervals of different encodings are never compared with each other, which is what
    distinguishes this rule from interval_direction applied to the pooled lists.
    """
    arms = sorted(set(intervals_by_arm_a) | set(intervals_by_arm_b))
    if not arms:
        return 0
    dirs = [interval_direction(intervals_by_arm_a.get(arm, []), intervals_by_arm_b.get(arm, []), eps)
            for arm in arms]
    if all(d == 1 for d in dirs):
        return 1
    if all(d == -1 for d in dirs):
        return -1
    return 0


def point_category(points_a, points_b, eps=EPS):
    """Classify the point comparison without collapsing ties into disagreement."""
    diffs = [a - b for a in points_a for b in points_b]
    if not diffs:
        return "empty", 0
    if all(d > eps for d in diffs):
        return "unanimous_positive", 1
    if all(d < -eps for d in diffs):
        return "unanimous_negative", -1
    if all(abs(d) <= eps for d in diffs):
        return "unanimous_tie", 0
    if any(d > eps for d in diffs) and any(d < -eps for d in diffs):
        return "conflicting", 0
    return "mixed_with_ties", 0


def agreement(predicted, measured):
    """None when either side does not resolve a direction; otherwise True or False."""
    if predicted == 0 or measured == 0:
        return None
    return predicted == measured
