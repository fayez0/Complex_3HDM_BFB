"""Exact positive-scale normalization before universal sign classification.

A coefficient vector p=scale*r, scale>0, has V_p=scale*V_r, so its quartic sign
is unchanged. The outer certificate retains BOTH p and the unmodified complete
decision for r. No tiny original-scale obstruction polynomials are expanded.
This addresses overall-scale growth, not every difficult relative hierarchy.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,time
from universal_regularizer import rational_data,KEYS
from robust_sign import classify as normalized_classify

def classify(data,max_halvings=256,max_shears=121):
    start=time.monotonic();p=rational_data(data);scale=max(abs(x) for k in KEYS for x in p[k]) or 1
    r={k:[str(x/scale) for x in p[k]] for k in KEYS};inner=normalized_classify(r,max_halvings,max_shears)
    return {'format':'positive-scale-normalized-v1','couplings':{k:list(map(str,p[k])) for k in KEYS},'positive_scale':str(scale),'normalized_decision':inner,'quartic_status':inner['quartic_status'],'original_deformation_radius':str(scale*__import__('sympy').Rational(inner['barrier']['sample_epsilon'])),'seconds':time.monotonic()-start,'scope':'Original sign decision via exact positive scaling; original singular minimum and minimizing fields not computed.'}

def verify_scale_link(record):
    if record.get('format')!='positive-scale-normalized-v1':raise ValueError('Unknown normalization record format')
    scale=Q(record['positive_scale'])
    if scale<=0:raise ValueError('Scale must be strictly positive')
    inner=record['normalized_decision']
    for k in KEYS:
        if len(record['couplings'][k])!=3 or len(inner['couplings'][k])!=3:raise ValueError('Wrong coefficient count')
        if any(Q(x)!=scale*Q(y) for x,y in zip(record['couplings'][k],inner['couplings'][k])):raise ValueError('Original/normalized coefficient mismatch')
    if record['quartic_status']!=inner['quartic_status']:raise ValueError('Original/normalized sign mismatch')
    eta=scale*Q(inner['barrier']['sample_epsilon'])
    if Q(record['original_deformation_radius'])!=eta:raise ValueError('Wrong original deformation radius')
    return {'scale_identity_verified':True,'positive_scale':str(scale),'original_deformation_radius':str(eta)}

def audit(record):
    from verify_complete import audit_decision
    link=verify_scale_link(record);inner=record['normalized_decision'];special=inner['barrier'].get('path_construction')=='anchor_homogeneity'
    if special:
        from homogeneous_anchor_path import audit as audit_path
        path_audit=audit_path(inner['barrier'])
    proof=audit_decision(inner,check_barrier=not special)
    if special:proof['barrier']=path_audit
    return {'verified':True,'scale_link':link,'normalized_complete_audit':proof,'status':record['quartic_status'],'scope':'Same completeness-proof dependence as the original auditor; exact positive-scale link checked separately.'}

def physical_certificate(record,check=True):
    verify_scale_link(record)
    if record['quartic_status']=='unstable':
        from physical_certificates import rational_witness
        d={'quartic_status':'unstable','couplings':record['couplings'],'samples':record['normalized_decision']['samples']};out=rational_witness(d);out['source']='Rational original fields from a negative normalized positive-sample candidate; final negativity checked directly on original coefficients';return out
    if record['quartic_status']=='strictly_positive':
        if check:audit(record)
        eta=Q(record['original_deformation_radius'])
        return {'couplings':record['couplings'],'quartic_status':'strictly_positive','normalized_lower_bound':str(5*eta/6),'componentwise_robustness_radius':str(5*eta/8),'complete_decision_audit_performed_here':check,'proof_dependency':'Positive scale link plus the complete negative-sample positivity audit. Bound is not the exact minimum.'}
    raise ValueError('No nonzero two-sided margin is asserted at a flat boundary')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('input');ap.add_argument('output');ap.add_argument('--audit',action='store_true');g=ap.add_mutually_exclusive_group();g.add_argument('--verify',action='store_true');g.add_argument('--physical',action='store_true');args=ap.parse_args();raw=json.loads(Path(args.input).read_text())
    if args.verify:out=audit(raw)
    elif args.physical:out=physical_certificate(raw)
    else:
        out=classify(raw)
        if args.audit:out['complete_audit']=audit(out)
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n');print(out.get('quartic_status',out.get('status',out.get('normalized_value','verified'))))
