import json
from pathlib import Path
import cobra
import numpy as np
from src.cpu_readout import select_flux

ROOT=Path(__file__).resolve().parents[1]
CONFIG=json.loads((ROOT/'configs/study.json').read_text())


def model():
    m=cobra.Model('diamond');m.solver='glpk'
    a=cobra.Metabolite('a_c');b=cobra.Metabolite('b_c');c=cobra.Metabolite('c_c')
    def add(rid, stoich, bounds):
        r=cobra.Reaction(rid);r.add_metabolites(stoich);r.bounds=bounds;m.add_reactions([r])
    add('IN',{a:1},(1,1));add('R1',{a:-1,b:1},(0,10));add('R2',{a:-1,c:1},(0,10))
    add('EX_b',{b:-1},(0,10));add('EX_c',{c:-1},(0,10))
    return m


def test_joint_lexicographic_solution_and_uncertainty():
    m=model();w={r.id:(1 if r.id in ['R1','R2'] else 0) for r in m.reactions}
    out=select_flux(m,w,['EX_b','EX_c'],CONFIG)
    assert out['ranges']['EX_b'][1]-out['ranges']['EX_b'][0]>.9
    assert abs(out['point']['EX_b']) <= CONFIG['coordinate_fix_allowance'] + 1e-12
    assert abs(out['point']['EX_c']-1) <= CONFIG['coordinate_fix_allowance'] + 1e-12
    assert out['audit']['max_mass_balance_residual']<1e-7
    out2=select_flux(model(),w,['EX_c','EX_b'],CONFIG)
    assert max(abs(out['point'][r]-out2['point'][r]) for r in out['point'])<1e-6


def test_flux_units_preserve_canonical_point():
    m=model();w={r.id:(1 if r.id in ['R1','R2'] else 0) for r in m.reactions}
    base=select_flux(m,w,['EX_b','EX_c'],CONFIG)
    m=model();r=m.reactions.get_by_id('EX_c');a=10
    r.add_metabolites({met:coeff*(1/a-1) for met,coeff in r.metabolites.items()});r.bounds=(r.lower_bound*a,r.upper_bound*a)
    sec={r.id:1 for r in m.reactions};sec['EX_c']=1/a
    out=select_flux(m,w,['EX_b','EX_c'],CONFIG,coordinate_factors={'EX_c':1/a},secondary_costs=sec)
    assert max(abs(base['point'][r]-out['point'][r]) for r in base['point'])<1e-6
