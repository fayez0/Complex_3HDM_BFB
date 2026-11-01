"""Recheck a zero-barrier certificate without importing its constructor/classifier.

Reconstructs stationary rows by direct angular-matrix multiplication, NOT by
the first-harmonic cancellation used by the constructor. Rebuilds all 35 path
polynomials, verifies factorizations and every exact Bernstein inequality.
This is an independent implementation path, not formal proof verification.
"""
from __future__ import annotations
import json,time
from pathlib import Path
import sympy as s
u,v,e=s.symbols('u v epsilon'); edges=((0,1),(1,2),(2,0)); opposite=(1,2,0)
keys=('a','b','c','d_real','d_imag')
anchor={'a':[5,6,7],'b':[0,0,0],'c':[2,3,4],'d_real':[1,2,-1],'d_imag':[1,-1,2]}


def parse(text): return s.sympify(text,locals={'u':u,'v':v,'epsilon':e})


def direct_rows(p,sector):
    a,b,c,dr,di=[p[k] for k in keys];D=(1+u*u)*(1+v*v)
    cos=[(1-u*u)*(1+v*v),(1-u*u)*(1-v*v)+4*u*v,(1-v*v)*(1+u*u)]
    sin=[-2*u*(1+v*v),2*u*(1-v*v)-2*v*(1-u*u),2*v*(1+u*u)]
    w=[s.expand(dr[j]*sin[j]+di[j]*cos[j]) for j in range(3)]
    q=[s.expand(dr[j]*cos[j]-di[j]*sin[j]) for j in range(3)]
    A=s.diag(*a)
    for j,(i,k) in enumerate(edges): A[i,k]=A[k,i]=b[j]
    h=[w[j] for j in opposite]
    if sector=='neutral':
        m=list(A*s.Matrix(h))
        for j,(i,k) in enumerate(edges):
            m[i]+=(c[j]*D+q[j])*h[k]/D; m[k]+=(c[j]*D+q[j])*h[i]/D
        return [s.cancel(z) for z in m]
    delta=[s.expand(c[j]**2-dr[j]**2-di[j]**2) for j in range(3)];J=s.prod(delta)
    H=[s.cancel(-h[i]*J/delta[opposite[i]]) for i in range(3)];K=[0,0,0]
    for j,(i,k) in enumerate(edges):
        K[i]+=(c[j]*D-q[j])*H[k]/(delta[j]*D);K[k]+=(c[j]*D-q[j])*H[i]/(delta[j]*D)
    return [s.cancel(z) for z in A*s.Matrix(K)+s.Matrix(H)]


def resultant_fixed22(f,g,x):
    # This alternative uses the actual Sylvester matrix with fixed row degrees.
    F,G=s.Poly(f,x),s.Poly(g,x);a,b,c=[F.nth(k) for k in (2,1,0)];d,h,j=[G.nth(k) for k in (2,1,0)]
    # Laplace expansion of the 4x4 determinant, polynomial in its entries.
    return s.expand(a*a*j*j-a*b*h*j-2*a*c*d*j+a*c*h*h+b*b*d*j-b*c*d*h+c*c*d*d)


def reconstructed_factors(p):
    a,b,c,dr,di=[p[k] for k in keys];out={};A=s.diag(*a);X=s.Symbol('X');reports={}
    out['cycle']=s.expand(dr[0]*dr[1]*di[2]+dr[0]*di[1]*dr[2]+di[0]*dr[1]*dr[2]-di[0]*di[1]*di[2])
    for j,(i,k) in enumerate(edges): A[i,k]=A[k,i]=b[j];out[f'delta_{j}']=s.expand(c[j]**2-dr[j]**2-di[j]**2)
    for sector in ('neutral','charged'):
        m=direct_rows(p,sector);f=s.expand(m[0]-m[2]);g=s.expand(m[1]-m[2])
        out[f'{sector}_resultant_v_at_u0']=resultant_fixed22(f.subs(u,0),g.subs(u,0),v)
        out[f'{sector}_resultant_u_at_v0']=resultant_fixed22(f.subs(v,0),g.subs(v,0),u)
        out[f'{sector}_u_infinity']=resultant_fixed22(s.Poly(f,u).nth(2),s.Poly(g,u).nth(2),v)
        out[f'{sector}_v_infinity']=resultant_fixed22(s.Poly(f,v).nth(2),s.Poly(g,v).nth(2),u)
        domain=s.QQ.poly_ring(e);div=s.Poly(1+u*u,u,domain=domain)
        R=s.Poly(resultant_fixed22(f,g,v),u,domain=domain).exquo(div);S=s.Poly(resultant_fixed22(f,m[0],v),u,domain=domain).exquo(div)
        if (R.degree(),S.degree())!=(6,6): raise AssertionError('Eliminant degree mismatch')
        cr,r=R.primitive();cs,z=S.primitive();out[f'zero_{sector}_full']=s.expand(cr**6*cs**6*r.resultant(z));reports[sector]=(R,S)
    for i in range(3):out[f'zero_singleton_{i}']=a[i]
    for j,(i,k) in enumerate(edges):
        rho2=dr[j]**2+di[j]**2;h=b[j]+c[j];B=A.copy();B[i,k]=B[k,i]=h+X
        out[f'orthogonal_pair_{j}']=a[i]+a[k]-2*b[j]
        out[f'neutral_pair_norm_{j}']=(a[i]+a[k]-2*h)**2-4*rho2
        for label,expression in [('split_norm',sum(B.adjugate())),('zero_split_norm',B.det())]:
            r=s.Poly(expression,X).rem(s.Poly(X*X-rho2,X));out[f'{label}_{j}']=s.expand(r.nth(0)**2-r.nth(1)**2*rho2)
        out[f'zero_orthogonal_pair_{j}']=a[i]*a[k]-b[j]*b[j]
        out[f'zero_neutral_pair_norm_{j}']=(a[i]*a[k]-h*h-rho2)**2-4*h*h*rho2
    return {k:s.Poly(x,e,domain=s.QQ)for k,x in out.items()},reports


def verify(certificate: dict, reconstruct: bool=True) -> dict:
    start=time.monotonic();c=certificate;radius=s.Rational(c['interval_certificate']['radius']);sample=s.Rational(c['sample_epsilon'])
    if not (0<sample<=radius):raise AssertionError('Sample outside certified interval')
    if c['anchor_certificate']['anchor']!=anchor:raise AssertionError('Unrecognized regularizer')
    p={k:[s.Rational(x)+e*y for x,y in zip(c['couplings'][k],anchor[k])]for k in keys}
    recorded={k:s.Poly(parse(x),e,domain=s.QQ)for k,x in c['path_polynomials'].items()}
    if len(recorded)!=35:raise AssertionError('Missing barriers')
    if reconstruct:
        rebuilt,reports=reconstructed_factors(p)
        if recorded!=rebuilt:raise AssertionError('Path-polynomial reconstruction failed')
        for sector,(R,S) in reports.items():
            if s.Poly(parse(c['eliminants'][sector]['R']),u)!=s.Poly(R.as_expr(),u) or s.Poly(parse(c['eliminants'][sector]['S']),u)!=s.Poly(S.as_expr(),u):raise AssertionError('Eliminant record mismatch')
    checkmap={row['polynomial']:row for row in c['interval_certificate']['factor_interval_checks']};used=set()
    for name,poly in recorded.items():
        dec=c['interval_certificate']['decompositions'][name];product=s.Poly(s.Rational(dec['constant'])*e**dec['zero_order'],e,domain=s.QQ)
        for item in dec['factors']:
            f=s.Poly(parse(item['polynomial']),e,domain=s.QQ);product*=f**item['multiplicity'];used.add(item['polynomial'])
            if item['polynomial'] not in checkmap:raise AssertionError('Missing factor interval certificate')
        if product!=poly:raise AssertionError('Factorization mismatch: '+name)
    inequalities=0
    for key in used:
        f=s.Poly(parse(key),e,domain=s.QQ);n=f.degree();constant=f.nth(0)
        if not constant:raise AssertionError('Unremoved zero at origin')
        sign=int(s.sign(constant))
        for side in (-1,1):
            row=checkmap[key]['checks'][str(side)];expected=[]
            for k in range(n+1):
                b=sum(f.nth(j)*(side*radius)**j*s.binomial(k,j)/s.binomial(n,j)for j in range(k+1));expected.append(str(b))
                if sign*b<=0:raise AssertionError('Bernstein sign failed')
                inequalities+=1
            if row['bernstein_coefficients']!=expected or row['constant_sign']!=sign:raise AssertionError('Stored Bernstein record mismatch')
    return {'verified':True,'polynomials_reconstructed':35 if reconstruct else 0,'factorizations_checked':35,'nonzero_factors_checked':len(used),'exact_Bernstein_inequalities':inequalities,'radius':str(radius),'seconds':time.monotonic()-start,'scope':'Independent construction and interval arithmetic path; not independent verification of the stationary-completeness theorem.'}


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('certificate');ap.add_argument('--output');args=ap.parse_args();raw=json.loads(Path(args.certificate).read_text());c=raw.get('barrier',raw);out=verify(c);text=json.dumps(out,indent=2);print(text)
    if args.output:Path(args.output).write_text(text+'\n')
