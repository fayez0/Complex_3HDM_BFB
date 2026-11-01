"""Standalone negative-quartic witness checker using ONLY Python's standard library.

The input is the original potential and a 2-by-3 matrix of complex rational
Higgs fields. No root solver, stationary equations, or completeness theorem is
used. Results apply to the quartic, not the full quadratic-plus-quartic model.
"""
from __future__ import annotations
from fractions import Fraction as Q
from pathlib import Path
import argparse, json
KEYS=('a','b','c','d_real','d_imag'); EDGES=((0,1),(1,2),(2,0))

def rational(x):
    if isinstance(x,(bool,float)): raise ValueError('Use integers or exact rational strings')
    if not isinstance(x,(str,int,Q)): raise ValueError('Unsupported rational encoding')
    return Q(x)

def cadd(a,b): return (a[0]+b[0],a[1]+b[1])
def cmul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cconj(a): return (a[0],-a[1])

def evaluate(couplings,real,imag):
    """Return exact V4, norm S, and Gram matrix for arbitrary nonzero fields."""
    if set(couplings)!=set(KEYS) or any(len(couplings[k])!=3 for k in KEYS): raise ValueError('Expected exactly five length-three arrays')
    if len(real)!=2 or len(imag)!=2 or any(len(row)!=3 for row in real+imag): raise ValueError('Fields must be 2-by-3 real and imaginary arrays')
    a,b,c,dr,di=[[rational(z) for z in couplings[k]] for k in KEYS]
    x=[[(rational(real[h][i]),rational(imag[h][i])) for i in range(3)] for h in range(2)]
    gram=[[cadd(cmul(cconj(x[0][i]),x[0][j]),cmul(cconj(x[1][i]),x[1][j])) for j in range(3)] for i in range(3)]
    r=[gram[i][i][0] for i in range(3)];S=sum(r,Q(0))
    if S<=0: raise ValueError('The zero field is not an instability witness')
    v=sum((a[i]*r[i]*r[i]/2 for i in range(3)),Q(0))
    for e,(i,j) in enumerate(EDGES):
        z=gram[i][j];zz=cmul(z,z);v+=b[e]*r[i]*r[j]+c[e]*(z[0]**2+z[1]**2)+dr[e]*zz[0]-di[e]*zz[1]
    return v,S,gram

def verify(record):
    v,S,gram=evaluate(record['couplings'],record['fields_real'],record['fields_imag']);normalized=v/(S*S)
    if v>=0: raise ValueError('The supplied fields do not prove instability')
    for name,actual in (('quartic_value',v),('field_norm',S),('normalized_value',normalized)):
        if name in record and rational(record[name])!=actual: raise ValueError('Wrong recorded '+name)
    # At a fixed trace-one physical Gram matrix, the 15-component infinity
    # perturbation norm changes V4 by at most (4/3) times that norm.
    radius=-3*normalized/4
    if 'componentwise_robustness_radius' in record and rational(record['componentwise_robustness_radius'])!=radius: raise ValueError('Wrong robustness radius')
    return {'verified':True,'original_potential_unstable':True,'quartic_value':str(v),'field_norm':str(S),'normalized_value':str(normalized),'componentwise_robustness_radius':str(radius),'dependencies':'Python standard library only','interpretation':'Negative feasible value, not a claim of the exact minimum; robustness uses strict componentwise radius.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('witness');ap.add_argument('--output');args=ap.parse_args()
    result=verify(json.loads(Path(args.witness).read_text()));text=json.dumps(result,indent=2)
    if args.output: Path(args.output).write_text(text+'\n')
    print(text)
