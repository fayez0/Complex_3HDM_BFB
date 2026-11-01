"""Exact factored barriers for p=lambda*q, avoiding redundant symbolic growth.

Uses the homogeneous degrees already proved for the same 35 obstruction
polynomials. The audit rebuilds the anchor values through direct angular
matrices, checks p=lambda*q exactly, then checks every path identity and the
unchanged Bernstein certificate. It does not guess a sign from a small float.
"""
from __future__ import annotations
from functools import lru_cache
import sympy as s
from zero_barrier import T,anchor_zero_certificate,rational_data,ANCHOR,KEYS

def anchor_multiple(data):
    p=rational_data(data);lam=p['a'][0]/s.Rational(ANCHOR['a'][0])
    return lam if all(p[k][i]==lam*ANCHOR[k][i] for k in KEYS for i in range(3)) else None

def path_factors(data):
    lam=anchor_multiple(data)
    if lam is None:raise ValueError('Input is not a scalar multiple of the fixed anchor')
    anchor=anchor_zero_certificate();factors={name:s.Poly(s.Rational(row['value'])*(T+lam)**row['homogeneous_degree'],T) for name,row in anchor['factors'].items()};reports={}
    for sector,degree in (('neutral',8),('charged',20)):
        row=anchor['full_sector_eliminants'][sector]
        reports[sector]={key:str((T+lam)**degree*s.sympify(row[key])) for key in ('R','S')};reports[sector].update(degrees=[6,6],homogeneous_projection_degree=degree)
    return factors,reports

def degrees():
    d={'cycle':3}
    for e in range(3):
        d[f'delta_{e}']=2;d[f'orthogonal_pair_{e}']=1;d[f'neutral_pair_norm_{e}']=2;d[f'split_norm_{e}']=4;d[f'zero_singleton_{e}']=1;d[f'zero_orthogonal_pair_{e}']=2;d[f'zero_neutral_pair_norm_{e}']=4;d[f'zero_split_norm_{e}']=6
    for sector,n in (('neutral',8),('charged',20)):
        for tag in ('resultant_v_at_u0','resultant_u_at_v0','u_infinity','v_infinity'):d[sector+'_'+tag]=n
    d.update(zero_neutral_full=96,zero_charged_full=240);return d

@lru_cache(maxsize=1)
def independently_rebuilt_anchor():
    from verify_barrier import reconstructed_factors
    p={k:list(map(s.Rational,v)) for k,v in ANCHOR.items()};factors,reports=reconstructed_factors(p)
    if set(factors)!=set(degrees()) or any(f.is_zero or f.degree()!=0 for f in factors.values()):raise AssertionError('Anchor certificate failed')
    return factors,reports

def audit(certificate):
    from verify_barrier import verify,parse,u
    lam=anchor_multiple(certificate['couplings'])
    if lam is None:raise AssertionError('False anchor-multiple claim')
    actual,reports=independently_rebuilt_anchor();degree=degrees()
    if set(certificate['path_polynomials'])!=set(actual):raise AssertionError('Wrong obstruction list')
    for name,f in actual.items():
        expected=s.Poly(f.nth(0)*(T+lam)**degree[name],T)
        if s.Poly(parse(certificate['path_polynomials'][name]),T)!=expected:raise AssertionError('Wrong homogeneous path factor: '+name)
    for sector,n in (('neutral',8),('charged',20)):
        for key,anchor_poly in zip(('R','S'),reports[sector]):
            if s.Poly(parse(certificate['eliminants'][sector][key]),u,T)!=s.Poly((T+lam)**n*anchor_poly.as_expr(),u,T):raise AssertionError('Wrong homogeneous projection')
    result=verify(certificate,reconstruct=False);result.update(anchor_multiple_verified=str(lam),anchor_values_independently_rebuilt=35,homogeneous_path_identities_checked=35,scope='Direct anchor reconstruction plus exact homogeneous path identities and Bernstein checks; universal completeness proof remains a separate obligation.')
    return result
