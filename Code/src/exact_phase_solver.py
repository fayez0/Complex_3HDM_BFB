"""Exact algebraic phase enumeration for checkpoint 02.

The algorithm exhausts the finite tangent chart when its fiber tests succeed.
Nonlinear surviving fibers or singular stationary fibers are reported as unsupported,
never treated as absent. Exact arithmetic uses SymPy QQ algebraic extensions.
"""
from __future__ import annotations
from pathlib import Path
import json
import sympy as sp
from sympy.polys.domains import QQ
U,V,T=sp.symbols('u v t')

def expr(text): return sp.sympify(text,locals={'u':U,'v':V})
def rat(q): return sp.Rational(q.numerator,q.denominator)
def interval_add(a,b): return a[0]+b[0],a[1]+b[1]
def interval_mul(a,b):
    p=[x*y for x in a for y in b]; return min(p),max(p)
def interval_poly(coeffs,box):
    ans=(sp.S.Zero,sp.S.Zero)
    for c in coeffs: ans=interval_add(interval_mul(ans,box),(rat(c),rat(c)))
    return ans

class Fiber:
    def __init__(self,h):
        self.h=sp.Poly(h,U,domain=QQ); self.K=QQ.alg_field_from_poly(self.h,alias='a'); self.a=self.scalar(U)
    def scalar(self,f): return self.K.new([QQ.from_sympy(c) for c in sp.Poly(f,U,domain=QQ).rem(self.h).all_coeffs()])
    def polynomial(self,f,var=V): return sp.Poly.from_dict({powers:self.scalar(c) for powers,c in sp.Poly(f,var).terms()},var,domain=self.K)
    def as_expr(self,a):
        c=a.to_list(); return sp.Add(*[rat(x)*U**(len(c)-i-1) for i,x in enumerate(c)])
    def poly_expr(self,f): return sp.Add(*[self.as_expr(a)*f.gen**(f.degree()-i) for i,a in enumerate(f.rep.to_list())])
    def evaluate_poly(self,f,v):
        ans=self.K.zero
        for (i,j),c in sp.Poly(f,U,V).terms(): ans+=self.K.convert(QQ.from_sympy(c))*self.a**i*v**j
        return ans
    def evaluate(self,f,v):
        n,d=sp.fraction(sp.cancel(f)); return self.evaluate_poly(n,v)/self.evaluate_poly(d,v)
    def signs(self,values,box):
        box=tuple(box)
        for digits in (20,40,80,160,320,640):
            out={}; intervals={}; success=True
            for k,a in values.items():
                if not a: out[k]=0;intervals[k]=(sp.S.Zero,sp.S.Zero);continue
                lo,hi=interval_poly(a.to_list(),box);intervals[k]=(lo,hi)
                if lo>0: out[k]=1
                elif hi<0: out[k]=-1
                else: success=False;break
            if success: return out,intervals,box
            box=self.h.refine_root(*box,eps=sp.Rational(1,10**digits))
        raise ArithmeticError('Exact signs not separated within the resource budget')

def remove_zero_factors(G,polynomials):
    removed=[]
    for name,p in polynomials:
        while G.degree()>0:
            d=G.gcd(p)
            if d.degree()==0: break
            removed.append((name,d.degree()));G=G.exquo(d)
    return G.monic(),removed

def endpoint_audit(N,D):
    """Certify missing tangent axes and corner; reject any unresolved endpoint."""
    out={};x,y=U,V
    for label,x,y in [('u_infinity',U,V),('v_infinity',V,U)]:
        n=sp.Poly(N,x);d=sp.Poly(D,x);n0,n1=n.nth(4),n.nth(3);d0,d1=d.nth(4),d.nth(3)
        a=sp.Poly(n1*d0-n0*d1,y);b=sp.Poly(sp.diff(n0,y)*d0-n0*sp.diff(d0,y),y);g=a.gcd(b)
        poles=sp.Poly(n0,y).gcd(sp.Poly(d0,y));out[label]={'critical_gcd':str(g.as_expr()),'pole_gcd':str(poles.as_expr())}
        if (g.degree()>0 and g.count_roots(-sp.oo,sp.oo)) or (poles.degree()>0 and poles.count_roots(-sp.oo,sp.oo)): raise NotImplementedError('Real endpoint fiber requires a separate solve: '+label)
    n=sp.Poly(N,U,V);d=sp.Poly(D,U,V);n00=n.coeff_monomial(U**4*V**4);d00=d.coeff_monomial(U**4*V**4)
    a=n.coeff_monomial(U**3*V**4)*d00-n00*d.coeff_monomial(U**3*V**4);b=n.coeff_monomial(U**4*V**3)*d00-n00*d.coeff_monomial(U**4*V**3)
    out['corner']={'N':str(n00),'D':str(d00),'derivative_numerators':[str(a),str(b)]}
    if (d00!=0 and a==b==0) or (d00==0 and n00==0): raise NotImplementedError('Corner needs an explicit stationary solve')
    return out

def enumerate_regular(record,geometry,couplings,sector):
    N,D=expr(record['mu_numerator']),expr(record['mu_denominator']);f,g=[expr(record[k]) for k in ('du_numerator_factored','dv_numerator_factored')]
    r=[expr(x) for x in geometry['r']];z=[expr(x) for x in geometry['z']];weights=r if sector=='neutral' else z
    # Divide only factors certified to have nonzero active support. This reduces
    # algebraic work; final saturation below independently applies all exclusions.
    for k in (1,): f=sp.exquo(f,sp.gcd(f,sp.fraction(weights[k])[0]),U,V)
    for k in (2,): g=sp.exquo(g,sp.gcd(g,sp.fraction(weights[k])[0]),U,V)
    f=sp.Poly(f,U,V).primitive()[1].as_expr();g=sp.Poly(g,U,V).primitive()[1].as_expr()
    if sp.Poly(sp.gcd(f,g),U,V).total_degree()>0: raise NotImplementedError('A common stationary component requires a separate treatment')
    R=sp.resultant(f,g,V)
    if R==0: raise NotImplementedError('Zero resultant is not a no-root certificate')
    factors=sp.factor_list(R)[1];rows=[];candidates=[]; poles=[]
    excluded=[('phase_denominator',D)]+[(f'r{i+1}',sp.fraction(x)[0]) for i,x in enumerate(r)]
    if sector=='charged': excluded += [(f'z{i+1}',sp.fraction(x)[0]) for i,x in enumerate(z)]+[('nu',sp.fraction(sp.cancel(sum(z)))[0])]
    aa,bb,cc,dr,di=[couplings[k] for k in ('a','b','c','d_real','d_imag')];edges=[(0,1),(1,2),(2,0)]
    for h,multiplicity in factors:
        hp=sp.Poly(h,U);boxes=hp.intervals(eps=sp.Rational(1,10**14));row={'h':str(h),'multiplicity_in_resultant':multiplicity,'real_projection_roots':len(boxes)}
        if not boxes: rows.append(row);continue
        F=Fiber(h);G=F.polynomial(f).gcd(F.polynomial(g)).sqf_part().monic();row['fiber_gcd']=str(F.poly_expr(G));bad=G.gcd(F.polynomial(D))
        if bad.degree()>0: poles.append({'h':str(h),'gcd':str(F.poly_expr(bad))})
        clean,removed=remove_zero_factors(G,[(name,F.polynomial(p)) for name,p in excluded]);row['removed']=removed;row['surviving_fiber']=str(F.poly_expr(clean));rows.append(row)
        if clean.degree()==0: continue
        if clean.degree()!=1: raise NotImplementedError('Nonlinear surviving fiber requires extension-field real isolation')
        ve=-clean.rep.to_list()[1]/clean.rep.to_list()[0]
        if F.evaluate_poly(f,ve) or F.evaluate_poly(g,ve): raise AssertionError('Exact substitution failed')
        values={f'r{i+1}':F.evaluate(x,ve) for i,x in enumerate(r)};values.update(v=ve,mu=F.evaluate(N/D,ve));
        if sector=='charged':
            zz=[F.evaluate(x,ve) for x in z];nu=sum(zz,F.K.zero);values.update({f'z{i+1}':x for i,x in enumerate(zz)});values['nu']=nu;values.update({f's{i+1}':x/nu for i,x in enumerate(zz)})
            trace2=sum((values[f'r{i+1}']**2 for i in range(3)),F.K.zero)
            for e,(i,j) in enumerate(edges):
                q=F.evaluate(expr(geometry['q'][e]),ve);delta=cc[e]**2-dr[e]**2-di[e]**2;trace2+=2*zz[i]*zz[j]*(cc[e]**2+dr[e]**2+di[e]**2-2*cc[e]*q)/(delta**2)
            values['psd_margin']=F.K.one-trace2
        for index,(box,m) in enumerate(boxes):
            signs,intervals,refined=F.signs(values,box);physical=all(signs[f'r{i+1}']>0 for i in range(3))
            if sector=='charged': physical=physical and all(signs[f's{i+1}']>0 for i in range(3)) and signs['psd_margin']>=0
            item={'sector':sector,'h':str(h),'projection_root_index':index,'u_isolation':list(map(str,refined)),'v_polynomial':str(F.as_expr(ve)),'field_values':{k:str(F.as_expr(x)) for k,x in values.items()},'signs':signs,'physical':physical,'display_approximations':{k:float((a+b)/2) for k,(a,b) in intervals.items()},'value_interval':[str(intervals['mu'][0]/2),str(intervals['mu'][1]/2)]};candidates.append(item)
    return {'sector':sector,'reduced_phase_equations':[str(f),str(g)],'resultant_degree':sp.degree(R,U),'factor_audit':rows,'pole_fibers':poles,'regular_candidates':candidates,'endpoints':endpoint_audit(N,D)}

def run_fixture(folder):
    folder=Path(folder);j=json.loads((folder/'rational_fixture_polynomials.json').read_text());geo=json.loads((folder/'geometry.json').read_text());out={}
    for name,rec in j['sectors'].items():
        geometry=dict(geo['sectors'][name],q=geo['q']);out[name]=enumerate_regular(rec,geometry,j['couplings'],name);print(name,'exact regular roots',len(out[name]['regular_candidates']),'physical',sum(x['physical'] for x in out[name]['regular_candidates']),flush=True)
        (folder/'regular_exact_certificate.json').write_text(json.dumps(out,indent=2,default=int)+'\n')
    return out
if __name__=='__main__':
    import sys
    run_fixture(sys.argv[1] if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'validation')
