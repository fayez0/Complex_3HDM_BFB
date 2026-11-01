"""All-rational-input quartic sign decision through a certified zero-free interval.

The mathematical coverage uses the proposed finite-current completeness proof.
This is not independently formally verified code, and explicit arithmetic budget
errors must not be interpreted as a physical classification.
"""
from __future__ import annotations
import json,time
from pathlib import Path
import sympy as s
from zero_barrier import build_certificate,rational_data,ANCHOR,KEYS
from exact_current import sparse
from separating_current import solve_sector


def finite_sample(data: dict, max_shears: int=121) -> dict:
    """A negative feasible sparse value suffices; otherwise exhaust both sectors."""
    out={'couplings':data,'sparse':sparse(data),'sectors':{}}
    bad=[row for row in out['sparse'] if row.get('physical') and row['value_sign']<0]
    if bad:
        out.update(quartic_status='unstable',complete_enumeration=False,negative_witness_label=bad[0]['label'],reason='A feasible negative sparse configuration suffices for instability'); return out
    for sector in ('neutral','charged'): out['sectors'][sector]=solve_sector(data,sector,max_shears)
    signs=[r['value_sign'] for r in out['sparse'] if r.get('physical')]+[r['signs']['mu'] for sec in out['sectors'].values() for r in sec['candidates'] if r['physical']]
    out.update(quartic_status='unstable' if min(signs)<0 else ('nonnegative_with_flat_direction' if min(signs)==0 else 'strictly_positive'),complete_enumeration=True)
    return out


def classify(data: dict, barrier: dict|None=None, max_shears: int=121) -> dict:
    start=time.monotonic(); p=rational_data(data)
    supplied_audit=None
    if barrier is None: barrier=build_certificate(data)
    else:
        from verify_barrier import verify
        supplied_audit=verify(barrier)
    if {k:list(map(str,v)) for k,v in p.items()}!=barrier['couplings']: raise ValueError('Barrier input mismatch')
    epsilon=s.Rational(barrier['sample_epsilon'])
    if epsilon<=0: raise ValueError('Sample epsilon must be positive')
    out={'couplings':barrier['couplings'],'scope':'Exact three-way quartic sign classification at rational input using a root-free deformation interval; original minimum and original attaining fields not computed','barrier':barrier,'samples':{}}
    for side in (1,-1):
        sample={k:[str(x+side*epsilon*y) for x,y in zip(p[k],ANCHOR[k])] for k in KEYS}; st=time.monotonic()
        result=finite_sample(sample,max_shears); result['seconds']=time.monotonic()-st; out['samples'][str(side)]=result
        if result['quartic_status']=='nonnegative_with_flat_direction': raise AssertionError('A zero sample contradicts the certified exclusion interval')
    pair=tuple(out['samples'][str(k)]['quartic_status'] for k in (1,-1))
    mapping={('strictly_positive','strictly_positive'):'strictly_positive',('strictly_positive','unstable'):'nonnegative_with_flat_direction',('unstable','unstable'):'unstable'}
    if pair not in mapping: raise AssertionError('Sample signs contradict positive-regularizer monotonicity')
    out['quartic_status']=mapping[pair];out['seconds']=time.monotonic()-start
    if supplied_audit is not None:out['supplied_barrier_verification']=supplied_audit
    out['limitations']=['Depends on the written proposed stationary-completeness proof; no independent specialist/formal verification claimed.','Exact original minimum and original attaining fields are not returned.','Resource-budget or algebraic-backend failures raise exceptions, never a guessed physical verdict.']
    return out


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('input');ap.add_argument('output');ap.add_argument('--barrier');args=ap.parse_args()
    data=json.loads(Path(args.input).read_text()); barrier=json.loads(Path(args.barrier).read_text()) if args.barrier else None
    out=classify(data,barrier);Path(args.output).write_text(json.dumps(out,indent=2)+'\n');print(out['quartic_status'],'epsilon',out['barrier']['sample_epsilon'],'seconds',out['seconds'])
