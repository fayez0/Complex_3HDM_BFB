"""Universal sign wrapper with a scale-independent, exact radius fallback.

The original radius selector and classifier are unchanged. If their declared
256-halving radius budget is exhausted, a coefficient-dominance bound supplies
a rational radius and the SAME Bernstein/audit certificate. No root-solving
failure or other arithmetic exception is converted into a physical verdict.
"""
from __future__ import annotations
from pathlib import Path
import argparse,json,time
import sympy as s
from zero_barrier import path_factors,certified_radius,strip_zero,no_root_test,anchor_zero_certificate,rational_data
from universal_sign import classify as classify_with_barrier,finite_sample
from homogeneous_anchor_path import anchor_multiple,path_factors as anchor_path_factors,audit as audit_anchor_path

def dominance_radius(factors):
    """For h=c0+sum cj*t^j, use min(1, |c0|/(2*sum_{j>0}|cj|))."""
    unique={};decomposition={}
    for name,p in factors.items():
        power,h=strip_zero(p);constant,parts=s.factor_list(h);row={'zero_order':power,'factors':[]}
        for f,multiplicity in parts:
            if f.degree()==0:continue
            constant*=f.LC()**multiplicity;f=f.monic();key=str(f.as_expr());unique[key]=f;row['factors'].append({'polynomial':key,'multiplicity':multiplicity})
        row['constant']=str(constant);decomposition[name]=row
    radius=s.S.One;bounds=[]
    for key,f in unique.items():
        c0=abs(f.nth(0));tail=sum(abs(f.nth(j)) for j in range(1,f.degree()+1))
        if c0==0:raise AssertionError('Zero root was not stripped')
        bound=min(s.S.One,c0/(2*tail)) if tail else s.S.One;radius=min(radius,bound)
        bounds.append({'polynomial':key,'constant_magnitude':str(c0),'coefficient_tail_sum':str(tail),'guaranteed_radius':str(bound)})
    checks=[]
    for key,f in unique.items():
        row={'polynomial':key,'checks':{}}
        for side in (-1,1):
            ok,proof=no_root_test(f,radius,side)
            if not ok:raise AssertionError('Coefficient-bound radius failed its independent Bernstein acceptance test')
            row['checks'][str(side)]=proof
        checks.append(row)
    return radius,{'method':'exact coefficient-dominance radius with the unchanged Bernstein acceptance certificate','radius':str(radius),'halvings':0,'distinct_nonconstant_factors':len(unique),'decompositions':decomposition,'factor_interval_checks':checks,'coefficient_bounds':bounds}

def build_barrier(data,max_halvings=256,force_dominance=False):
    start=time.monotonic();lam=anchor_multiple(data);polys,eliminants=(path_factors(data) if lam is None else anchor_path_factors(data));fallback=None
    if force_dominance:radius,proof=dominance_radius(polys);fallback='explicitly requested coefficient bound'
    else:
        try:radius,proof=certified_radius(polys,max_halvings)
        except ArithmeticError as exc:
            if str(exc)!='Root-free interval budget exceeded':raise
            fallback=str(exc);radius,proof=dominance_radius(polys)
    out={'couplings':{k:list(map(str,v)) for k,v in rational_data(data).items()},'scope':'Exact punctured zero-energy exclusion interval; no sample classification yet','anchor_certificate':anchor_zero_certificate(),'path_polynomials':{k:str(v.as_expr()) for k,v in polys.items()},'eliminants':eliminants,'interval_certificate':proof,'sample_epsilon':str(radius),'seconds':time.monotonic()-start,'radius_fallback_reason':fallback}
    if lam is not None:out['path_construction']='anchor_homogeneity';out['anchor_multiple']=str(lam)
    return out

def classify(data,max_halvings=256,max_shears=121,force_dominance=False):
    barrier=build_barrier(data,max_halvings,force_dominance)
    if barrier.get('path_construction')!='anchor_homogeneity':return classify_with_barrier(data,barrier=barrier,max_shears=max_shears)
    start=time.monotonic();verification=audit_anchor_path(barrier);p=rational_data(data);eta=s.Rational(barrier['sample_epsilon']);samples={}
    from universal_regularizer import ANCHOR,KEYS
    for side in (1,-1):
        sample={k:[str(x+side*eta*y) for x,y in zip(p[k],ANCHOR[k])] for k in KEYS};samples[str(side)]=finite_sample(sample,max_shears)
    pair=tuple(samples[str(side)]['quartic_status'] for side in (1,-1));table={('strictly_positive','strictly_positive'):'strictly_positive',('strictly_positive','unstable'):'nonnegative_with_flat_direction',('unstable','unstable'):'unstable'}
    if pair not in table:raise AssertionError('Samples contradict the exact zero-free interval')
    return {'couplings':barrier['couplings'],'scope':'Same finite-sample sign criterion with an exact homogeneous anchor-path construction','barrier':barrier,'samples':samples,'quartic_status':table[pair],'supplied_barrier_verification':verification,'seconds':time.monotonic()-start}

def audit(record):
    from verify_complete import audit_decision
    special=record['barrier'].get('path_construction')=='anchor_homogeneity'
    if special:path_audit=audit_anchor_path(record['barrier'])
    result=audit_decision(record,check_barrier=not special)
    if special:result['barrier']=path_audit
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('input');ap.add_argument('output');ap.add_argument('--force-dominance',action='store_true');args=ap.parse_args()
    out=classify(json.loads(Path(args.input).read_text()),force_dominance=args.force_dominance);Path(args.output).write_text(json.dumps(out,indent=2)+'\n');print(out['quartic_status'],out['barrier']['sample_epsilon'])
