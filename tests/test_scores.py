import numpy as np
import cobra
from scipy.stats import rankdata
from src.cpu_model import score_weights


def fixture():
    m=cobra.Model("score_fixture")
    a,b=cobra.Metabolite("a_c"),cobra.Metabolite("b_c")
    for i in range(6):
        r=cobra.Reaction(f"r{i}");r.add_metabolites({a:-1,b:1});m.add_reactions([r])
    r=cobra.Reaction("EX_b");r.add_metabolites({b:-1});m.add_reactions([r])
    return m


def test_monotone_grid_preserves_true_ties_and_support():
    m=fixture();ids=[f"r{i}" for i in range(5)];a=np.array([0.,1.,1.,2.,100.]);k=3.
    z=a/(a+k);base_ranks=rankdata(z,method="average")
    costs=[]
    for gamma in [.5,.8,1.,1.25,2.]:
        assert np.array_equal(rankdata(z**gamma,method="average"),base_ranks)
        w=score_weights(m,ids,a,k,f"magnitude_g{gamma:.2f}",.001)
        assert w["r1"]==w["r2"]
        assert w["r5"]==1. and w["EX_b"]==0.
        assert np.isclose(np.mean([w[r] for r in ids]),1.)
        costs.append([w[r] for r in ids])
    assert np.max(np.abs(np.array(costs)[0]-np.array(costs)[-1]))>.1
    o=score_weights(m,ids,a,k,"ordinal",.001)
    assert o["r1"]==o["r2"]


def test_expression_and_scale_unit_change_preserves_costs():
    m=fixture();ids=[f"r{i}" for i in range(5)];a=np.array([0.,1.,1.,2.,100.])
    for arm in ["ordinal","magnitude_g1.00","magnitude_g0.50","uniform_pfba"]:
        w=score_weights(m,ids,a,3.,arm,.001)
        v=score_weights(m,ids,a*1000.,3000.,arm,.001)
        assert max(abs(w[r]-v[r]) for r in w)<1e-14
