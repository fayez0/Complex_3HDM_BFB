"""Exact fixed-direction regularization for the complex Weinberg quartic.

This module constructs homogeneous algebraic obstructions and their anchor
certificate. It does NOT implement signs of arbitrary Puiseux-algebraic roots.
All original checkpoint-02 modules are preserved unchanged.
"""
from __future__ import annotations
import json
from pathlib import Path
import sympy as s
U,V=s.symbols('u v'); E=((0,1),(1,2),(2,0)); OPP=(1,2,0)
ANCHOR={'a':[5,6,7],'b':[0,0,0],'c':[2,3,4],'d_real':[1,2,-1],'d_imag':[1,-1,2]}
KEYS=('a','b','c','d_real','d_imag')

def rational_data(data):
    if any(len(data[k])!=3 for k in KEYS): raise ValueError('Each coefficient array must have length three')
    if any(isinstance(v,float) for k in KEYS for v in data[k]): raise ValueError('Use integers or rational strings, not binary floats')
    return {k:[s.Rational(v) for v in data[k]] for k in KEYS}

def polynomial_geometry(data,sector):
    """Unscaled homogeneous phase equations; degree 2 (N), 5 (C) in couplings."""
    a,b,c,dr,di=[list(map(s.sympify,data[k])) for k in KEYS]
    co=[(1-U**2)*(1+V**2),(1-U**2)*(1-V**2)+4*U*V,(1-V**2)*(1+U**2)]
    si=[-2*U*(1+V**2),2*U*(1-V**2)-2*V*(1-U**2),2*V*(1+U**2)]
    w=[dr[e]*si[e]+di[e]*co[e] for e in range(3)]; h=s.Matrix([w[e] for e in OPP]); extra=[]
    for opp in OPP:
        j,k=[e for e in range(3) if e!=opp]
        extra.append((dr[j]*di[k]+di[j]*dr[k])*co[opp]-(dr[j]*dr[k]-di[j]*di[k])*si[opp])
    A=s.diag(*a)
    for e,(i,j) in enumerate(E): A[i,j]=A[j,i]=b[e]
    delta=[c[e]**2-dr[e]**2-di[e]**2 for e in range(3)]
    if sector=='neutral':
        B0=A.copy()
        for e,(i,j) in enumerate(E): B0[i,j]=B0[j,i]=b[e]+c[e]
        m=B0*h+s.Matrix(extra)
    elif sector=='charged':
        H=[];K=[]
        for i,opp in enumerate(OPP):
            j,k=[e for e in range(3) if e!=opp]
            H.append(-w[opp]*delta[j]*delta[k]);K.append(delta[opp]*(extra[i]-c[j]*w[k]-c[k]*w[j]))
        m=A*s.Matrix(K)+s.Matrix(H)
    else: raise ValueError('sector must be neutral or charged')
    return s.expand(m[0]-m[2]),s.expand(m[1]-m[2])

def resultant22(f,g,x):
    """Fixed-degree Sylvester determinant, including degree-drop factors."""
    F,G=s.Poly(f,x),s.Poly(g,x)
    if max(F.degree(),G.degree())>2: raise ValueError('Expected quadratics or lower')
    a,b,c=[F.nth(k) for k in (2,1,0)];d,e,h=[G.nth(k) for k in (2,1,0)]
    return s.expand((a*h-c*d)**2-(a*e-b*d)*(b*h-c*e))

def obstructions(data):
    a,b,c,dr,di=[list(map(s.sympify,data[k])) for k in KEYS]; out={}
    out['cycle']=(s.expand(dr[0]*dr[1]*di[2]+dr[0]*di[1]*dr[2]+di[0]*dr[1]*dr[2]-di[0]*di[1]*di[2]),3)
    for e in range(3): out[f'delta_{e}']=(s.expand(c[e]**2-dr[e]**2-di[e]**2),2)
    for sector,degree in [('neutral',8),('charged',20)]:
        f,g=polynomial_geometry(data,sector)
        out[f'{sector}_resultant_v_at_u0']=(resultant22(f.subs(U,0),g.subs(U,0),V),degree)
        out[f'{sector}_resultant_u_at_v0']=(resultant22(f.subs(V,0),g.subs(V,0),U),degree)
        out[f'{sector}_u_infinity']=(resultant22(s.Poly(f,U).nth(2),s.Poly(g,U).nth(2),V),degree)
        out[f'{sector}_v_infinity']=(resultant22(s.Poly(f,V).nth(2),s.Poly(g,V).nth(2),U),degree)
    A=s.diag(*a);X=s.Symbol('X',real=True)
    for e,(i,j) in enumerate(E): A[i,j]=A[j,i]=b[e]
    for e,(i,j) in enumerate(E):
        rho2=dr[e]**2+di[e]**2
        out[f'orthogonal_pair_{e}']=(s.expand(a[i]+a[j]-2*b[e]),1)
        out[f'neutral_pair_norm_{e}']=(s.expand((a[i]+a[j]-2*(b[e]+c[e]))**2-4*rho2),2)
        B=A.copy();B[i,j]=B[j,i]=b[e]+c[e]+X;J=s.expand((s.ones(1,3)*B.adjugate()*s.ones(3,1))[0]);J=s.rem(J,X**2-rho2,X)
        j0,j1=s.Poly(J,X).nth(0),s.Poly(J,X).nth(1)
        out[f'split_norm_{e}']=(s.expand(j0*j0-j1*j1*rho2),4)
    return out

def anchor_certificate():
    p=rational_data(ANCHOR);obs=obstructions(p)
    assert all(value!=0 for value,_ in obs.values())
    # Uniform inequalities 5/6 <= V_anchor <= 7/2 on trace-one PSD matrices.
    for e in range(3): assert p['c'][e]>0 and p['c'][e]**2>p['d_real'][e]**2+p['d_imag'][e]**2
    assert min(p['a'])==5 and max(p['a'])==7 and p['b']==[0]*3
    assert max(p['c'][e]+abs(p['d_real'][e])+abs(p['d_imag'][e]) for e in range(3))<=7
    out={'anchor':ANCHOR,'lower_bound':'5/6','upper_bound':'7/2','obstructions':{k:{'value':str(v),'homogeneous_degree':d} for k,(v,d) in obs.items()},'product_degree':sum(d for _,d in obs.values()),'product_value':str(s.prod(v for v,_ in obs.values()))}
    out['all_nonzero']=True;return out

if __name__=='__main__':
    target=Path(__file__).resolve().parents[1]/'validation'/'anchor_certificate.json';target.write_text(json.dumps(anchor_certificate(),indent=2)+'\n');print(target)
