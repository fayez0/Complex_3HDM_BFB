"""Exact zero-energy barriers for the complex Weinberg quartic.

Adds 14 homogeneous obstructions to the 21 checkpoint-03 regularity factors.
Their nonvanishing excludes ALL stationary configurations of energy zero in
that regular sector. The fixed positive anchor makes every path polynomial
nonzero. A punctured root-free interval then reduces the all-coupling SIGN
problem to exact finite rational samples; it does not compute the original
minimum or original attaining fields.

No floating point is used for polynomial construction or interval certificates.
Original checkpoint-03 modules are imported without modification.
"""
from __future__ import annotations
import json, time
from functools import lru_cache
from pathlib import Path
import sympy as s
from universal_regularizer import ANCHOR, KEYS, E, OPP, U, V, rational_data, obstructions, resultant22
T=s.Symbol('epsilon')


def stationary_rows(data: dict, sector: str) -> list:
    """Denominator-cleared homogeneous stationary rows, degrees 2 and 5."""
    a,b,c,dr,di=[list(map(s.sympify,data[k])) for k in KEYS]
    co=[(1-U**2)*(1+V**2),(1-U**2)*(1-V**2)+4*U*V,(1-V**2)*(1+U**2)]
    si=[-2*U*(1+V**2),2*U*(1-V**2)-2*V*(1-U**2),2*V*(1+U**2)]
    w=[dr[e]*si[e]+di[e]*co[e] for e in range(3)]; h=s.Matrix([w[e] for e in OPP]); extra=[]
    for opp in OPP:
        j,k=[e for e in range(3) if e!=opp]
        extra.append((dr[j]*di[k]+di[j]*dr[k])*co[opp]-(dr[j]*dr[k]-di[j]*di[k])*si[opp])
    A=s.diag(*a)
    for e,(i,j) in enumerate(E): A[i,j]=A[j,i]=b[e]
    if sector=='neutral':
        B=A.copy()
        for e,(i,j) in enumerate(E): B[i,j]=B[j,i]=b[e]+c[e]
        m=B*h+s.Matrix(extra)
    elif sector=='charged':
        delta=[c[e]**2-dr[e]**2-di[e]**2 for e in range(3)]; H=[]; K=[]
        for i,opp in enumerate(OPP):
            j,k=[e for e in range(3) if e!=opp]
            H.append(-w[opp]*delta[j]*delta[k]); K.append(delta[opp]*(extra[i]-c[j]*w[k]-c[k]*w[j]))
        m=A*s.Matrix(K)+s.Matrix(H)
    else: raise ValueError('Expected neutral or charged sector')
    return [s.expand(x) for x in m]


def full_zero_factor(data: dict, sector: str, parameter=None) -> tuple:
    """Return Z=Res_u(R,S), after removing ONLY the no-real-root factor u^2+1.

    R=Res_v(F,G)/(u^2+1), S=Res_v(F,m1)/(u^2+1).
    Fixed degrees are 2 internally and 6 externally. Path and anchor leading
    degrees are checked. A vanished result is an error, never a certificate.
    """
    m=stationary_rows(data,sector); f=s.expand(m[0]-m[2]); g=s.expand(m[1]-m[2])
    domain=s.QQ if parameter is None else s.QQ.poly_ring(parameter)
    div=s.Poly(1+U**2,U,domain=domain)
    R=s.Poly(resultant22(f,g,V),U,domain=domain).exquo(div)
    S=s.Poly(resultant22(f,m[0],V),U,domain=domain).exquo(div)
    if R.degree()!=6 or S.degree()!=6: raise ArithmeticError('Unexpected projection degree; fixed-degree proof needs audit')
    cr,r=R.primitive(); cs,z=S.primitive()
    value=s.expand(cr**6*cs**6*r.resultant(z))
    if value==0: raise ArithmeticError('Identically vanishing zero-energy obstruction')
    return value,{'R':str(R.as_expr()),'S':str(S.as_expr()),'R_content':str(cr),'S_content':str(cs),'degrees':[6,6]}


def sparse_zero_factors(data: dict) -> dict:
    a,b,c,dr,di=[list(map(s.sympify,data[k])) for k in KEYS]; A=s.diag(*a); X=s.Symbol('X'); out={}
    for e,(i,j) in enumerate(E): A[i,j]=A[j,i]=b[e]
    for i in range(3): out[f'zero_singleton_{i}']=(a[i],1)
    for e,(i,j) in enumerate(E):
        rho2=dr[e]**2+di[e]**2; h=b[e]+c[e]
        out[f'zero_orthogonal_pair_{e}']=(s.expand(a[i]*a[j]-b[e]**2),2)
        out[f'zero_neutral_pair_norm_{e}']=(s.expand((a[i]*a[j]-h*h-rho2)**2-4*h*h*rho2),4)
        B=A.copy(); B[i,j]=B[j,i]=h+X
        D=s.Poly(B.det(),X).rem(s.Poly(X**2-rho2,X)); d0,d1=D.nth(0),D.nth(1)
        out[f'zero_split_norm_{e}']=(s.expand(d0*d0-d1*d1*rho2),6)
    return out


@lru_cache(maxsize=1)
def anchor_zero_certificate() -> dict:
    q=rational_data(ANCHOR); factors=obstructions(q); factors.update(sparse_zero_factors(q)); reports={}
    for sector,degree in [('neutral',96),('charged',240)]:
        val,rep=full_zero_factor(q,sector); factors[f'zero_{sector}_full']=(val,degree); reports[sector]=rep
    if len(factors)!=35 or any(v==0 for v,d in factors.values()): raise AssertionError('Anchor obstruction check failed')
    return {'anchor':ANCHOR,'factor_count':len(factors),'product_degree':sum(d for v,d in factors.values()),'factors':{k:{'value':str(v),'homogeneous_degree':d} for k,(v,d) in factors.items()},'full_sector_eliminants':reports,'all_nonzero':True}


def path_factors(data: dict) -> tuple:
    """Construct all 35 nonzero univariate factors on p+epsilon*q."""
    data=rational_data(data); path={k:[x+T*y for x,y in zip(data[k],ANCHOR[k])] for k in KEYS}; factors=obstructions(path); factors.update(sparse_zero_factors(path)); reports={}
    for sector,degree in [('neutral',96),('charged',240)]:
        val,rep=full_zero_factor(path,sector,T); factors[f'zero_{sector}_full']=(val,degree); reports[sector]=rep
    anchor=anchor_zero_certificate(); out={}
    for name,(value,degree) in factors.items():
        p=s.Poly(value,T,domain=s.QQ)
        if p.is_zero or p.degree()!=degree or p.LC()!=s.Rational(anchor['factors'][name]['value']): raise AssertionError('Homogeneous anchor identity failed: '+name)
        out[name]=p
    return out,reports


def strip_zero(p: s.Poly) -> tuple:
    """Remove the exact zero-root power; never removes any nonzero root."""
    if p.is_zero: raise ValueError('Zero polynomial cannot certify an interval')
    n=min(k[0] for k in p.monoms()); return n,s.Poly.from_dict({(k[0]-n,):v for k,v in p.terms()},T,domain=s.QQ)


def bernstein_coefficients(p: s.Poly, radius, side: int) -> list:
    """Bernstein coefficients on epsilon in [0,r] or [-r,0], exact QQ."""
    n=p.degree(); coeff=[p.nth(j)*(side*radius)**j for j in range(n+1)]
    return [sum(coeff[j]*s.Rational(s.binomial(k,j),s.binomial(n,j)) for j in range(k+1)) for k in range(n+1)]


def no_root_test(p: s.Poly, radius, side: int) -> tuple:
    power,h=strip_zero(p); coeff=bernstein_coefficients(h,radius,side); sign=int(s.sign(h.nth(0)))
    return all(sign*x>0 for x in coeff),{'zero_order':power,'degree_after_zero_removal':h.degree(),'constant_sign':sign,'bernstein_coefficients':list(map(str,coeff))}


def certified_radius(factors: dict, max_halvings: int=256) -> tuple:
    """Factor each polynomial, deduplicate, and certify a common two-sided radius.

    Strict same-sign Bernstein coefficients prove nonvanishing on the entire
    CLOSED interval after its zero-root factor is removed. Halving must
    eventually succeed because every remaining polynomial is nonzero at zero.
    A resource-budget failure is explicit, not a BFB verdict.
    """
    unique={}; decomposition={}
    for name,p in factors.items():
        power,h=strip_zero(p); content,parts=s.factor_list(h); decomposition[name]={'zero_order':power,'constant':str(content),'factors':[]}
        for f,multiplicity in parts:
            if f.degree()==0: continue
            content*=f.LC()**multiplicity; f=f.monic(); key=str(f.as_expr()); unique[key]=f
            decomposition[name]['factors'].append({'polynomial':key,'multiplicity':multiplicity})
        decomposition[name]['constant']=str(content)
    radius=s.S.One; steps=0
    for f in unique.values():
        while not all(no_root_test(f,radius,side)[0] for side in (-1,1)):
            radius/=2; steps+=1
            if steps>max_halvings: raise ArithmeticError('Root-free interval budget exceeded')
    checks=[]
    for key,p in unique.items():
        row={'polynomial':key,'checks':{}}
        for side in (-1,1):
            ok,report=no_root_test(p,radius,side)
            if not ok: raise AssertionError('Lost a previously certified root-free interval')
            row['checks'][str(side)]=report
        checks.append(row)
    return radius,{'method':'exact Bernstein coefficient signs after removing epsilon powers','radius':str(radius),'halvings':steps,'distinct_nonconstant_factors':len(unique),'decompositions':decomposition,'factor_interval_checks':checks}


def build_certificate(data: dict, max_halvings: int=256) -> dict:
    start=time.monotonic(); polys,eliminants=path_factors(data); radius,proof=certified_radius(polys,max_halvings)
    return {'couplings':{k:list(map(str,v)) for k,v in rational_data(data).items()},'scope':'Exact punctured zero-energy exclusion interval; no sample classification yet','anchor_certificate':anchor_zero_certificate(),'path_polynomials':{k:str(v.as_expr()) for k,v in polys.items()},'eliminants':eliminants,'interval_certificate':proof,'sample_epsilon':str(radius),'seconds':time.monotonic()-start}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('input'); parser.add_argument('output'); args=parser.parse_args()
    data=json.loads(Path(args.input).read_text()); out=build_certificate(data); Path(args.output).write_text(json.dumps(out,indent=2)+'\n'); print('certified radius',out['sample_epsilon'])
