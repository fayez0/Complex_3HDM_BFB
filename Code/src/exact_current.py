"""Exact finite-current classifier prototype, checkpoint 02.

Supports rational inputs with nonreal cycle and nonsingular charged angular
blocks when the detected phase systems have finite, linear lifting fibers.
Unsupported cases raise explicit exceptions; no heuristic BFB answer is returned.
Sparse singular systems are covered by the boundary-descent lemma in the notes.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import sympy as s
from exact_phase_solver import U,V,Fiber,remove_zero_factors,interval_poly
E=((0,1),(1,2),(2,0));OPPOSITE=(1,2,0)

def qinput(value):
    if isinstance(value,float): raise ValueError('Use integers or exact rational strings, not binary floats')
    z=s.Rational(value)
    return z

def geometry(data,sector,chart='finite'):
    a,b,c,dr,di=[[qinput(v) for v in data[k]] for k in ('a','b','c','d_real','d_imag')]
    if any(len(z)!=3 for z in (a,b,c,dr,di)): raise ValueError('All input arrays must have length three')
    cycle=s.expand(s.prod(dr[e]+s.I*di[e] for e in range(3)))
    if s.im(cycle)==0: raise NotImplementedError('Real invariant-cycle stratum: rephase and use a real-coefficient solver')
    delta=[c[e]**2-dr[e]**2-di[e]**2 for e in range(3)]
    if sector=='charged' and any(x==0 for x in delta): raise NotImplementedError('Singular angular-block stratum is not implemented')
    C=lambda t:(1-t*t)/(1+t*t);S=lambda t:2*t/(1+t*t)
    if chart=='finite': cx,sx,cy,sy=C(U),S(U),C(V),S(V);den=(1+U**2)*(1+V**2)
    elif chart=='u_infinity': cx,sx,cy,sy=-1,0,C(U),S(U);den=1+U**2
    elif chart=='v_infinity': cx,sx,cy,sy=C(U),S(U),-1,0;den=1+U**2
    elif chart=='corner': cx,sx,cy,sy=-1,0,-1,0;den=s.S.One
    else: raise ValueError('Unknown phase chart')
    co=[cx,cx*cy+sx*sy,cy];si=[-sx,sx*cy-cx*sy,sy];q=[s.cancel(dr[e]*co[e]-di[e]*si[e]) for e in range(3)]
    ar,ai=(dr,di) if sector=='neutral' else ([-dr[e]/delta[e] for e in range(3)],[-di[e]/delta[e] for e in range(3)])
    w=[s.cancel(ar[e]*si[e]+ai[e]*co[e]) for e in range(3)];h=s.Matrix([w[e] for e in OPPOSITE]);extra=[]
    for i,opp in enumerate(OPPOSITE):
        others=[e for e in range(3) if e!=opp];j,k=others;pr=ar[j]*ar[k]-ai[j]*ai[k];im=ar[j]*ai[k]+ai[j]*ar[k];extra.append(s.cancel(im*co[opp]-pr*si[opp]))
    A=s.diag(*a)
    for e,(i,j) in enumerate(E): A[i,j]=A[j,i]=b[e]
    if sector=='neutral':
        C0=A.copy()
        for e,(i,j) in enumerate(E):C0[i,j]=C0[j,i]=b[e]+c[e]
        m=(C0*h+s.Matrix(extra)).applyfunc(s.cancel);norm=s.cancel(sum(h));r=(h/norm).applyfunc(s.cancel);z=s.zeros(3,1)
    else:
        C0=s.zeros(3)
        for e,(i,j) in enumerate(E):C0[i,j]=C0[j,i]=c[e]/delta[e]
        k=(C0*h+s.Matrix(extra)).applyfunc(s.cancel);m=(A*k+h).applyfunc(s.cancel);norm=s.cancel(sum(k));r=(k/norm).applyfunc(s.cancel);z=(-h/norm).applyfunc(s.cancel)
    f=s.Poly(s.cancel(den*(m[0]-m[2])),U,V).primitive()[1].as_expr();g=s.Poly(s.cancel(den*(m[1]-m[2])),U,V).primitive()[1].as_expr()
    return dict(f=f,g=g,r=list(r),z=list(z),mu=s.cancel(m[0]/norm),h=list(h),norm=norm,W=s.cancel(sum(h)),q=q,delta=delta,c=c,dr=dr,di=di)

def evaluate_candidates(h,F,G,geo,sector,chart):
    if G.degree()!=1: raise NotImplementedError('Nonlinear lifting fiber; a separating projection or extension-field solver is required')
    ve=-G.rep.to_list()[1]/G.rep.to_list()[0];v=ve if chart=='finite' else F.K.zero
    if F.evaluate_poly(geo['f'],v) or F.evaluate_poly(geo['g'],v): raise AssertionError('Exact lifted phase substitution failed')
    values={f'r{i+1}':F.evaluate(x,v) for i,x in enumerate(geo['r'])};values.update(mu=F.evaluate(geo['mu'],v))
    if sum((values[f'r{i+1}']for i in range(3)),F.K.zero)!=F.K.one: raise AssertionError('Exact normalization failed')
    if chart=='finite':values['v']=v
    if sector=='charged':
        zz=[F.evaluate(x,v) for x in geo['z']];nu=sum(zz,F.K.zero);values['nu']=nu;values.update({f's{i+1}':x/nu for i,x in enumerate(zz)});trace2=sum((values[f'r{i+1}']**2 for i in range(3)),F.K.zero)
        for e,(i,j) in enumerate(E):
            q=F.evaluate(geo['q'][e],v);c,dr,di,delta=[geo[k][e] for k in ('c','dr','di','delta')];trace2+=2*zz[i]*zz[j]*(c*c+dr*dr+di*di-2*c*q)/(delta*delta)
        values['psd_margin']=F.K.one-trace2
    candidates=[]
    for index,(box,mult) in enumerate(s.Poly(h,U).intervals(eps=s.Rational(1,10**40))):
        signs,ivals,box=F.signs(values,box);physical=all(signs[f'r{i+1}']>0 for i in range(3))
        if sector=='charged':physical=physical and all(signs[f's{i+1}']>0 for i in range(3)) and signs['psd_margin']>=0
        if chart=='finite':xy=[float(sum(box)/2),float(sum(ivals['v'])/2)]
        elif chart=='u_infinity':xy=['infinity',float(sum(box)/2)]
        elif chart=='corner':xy=['infinity','infinity']
        else:xy=[float(sum(box)/2),'infinity']
        candidates.append({'sector':sector,'chart':chart,'h':str(h),'root_index':index,'u_interval':list(map(str,box)),'v_polynomial':str(F.as_expr(v)),'field_values':{k:str(F.as_expr(x)) for k,x in values.items()},'signs':signs,'physical':bool(physical),'value_interval':[str(ivals['mu'][0]/2),str(ivals['mu'][1]/2)],'display_value':float(sum(ivals['mu'])/4),'display_tangents':xy,'display_quantities':{k:float(sum(z)/2) for k,z in ivals.items()}})
    return candidates

def solve_sector(data,sector):
    out={'sector':sector,'charts':{},'candidates':[]}
    for chart in ('finite','u_infinity','v_infinity','corner'):
        geo=geometry(data,sector,chart);f,g=geo['f'],geo['g'];report={'equations':[str(f),str(g)]}
        if chart=='corner':
            report['stationary']=f==g==0
            if report['stationary']:
                if geo['norm']==0 or geo['W']==0 or any(x==0 for x in geo['h']):report['excluded']='zero normalization or current'
                else:
                    F=Fiber(U);G=s.Poly.from_dict({(1,):F.K.one},V,domain=F.K);out['candidates']+=evaluate_candidates(U,F,G,geo,sector,chart)
            out['charts'][chart]=report;continue
        if chart=='finite':
            if s.Poly(s.gcd(f,g),U,V).total_degree()>0:raise NotImplementedError('A positive-dimensional stationary set requires component treatment')
            R=s.resultant(f,g,V)
            if R==0:raise NotImplementedError('Vanishing resultant requires a different projection')
            factors=s.factor_list(R)[1];report['resultant']=str(s.factor(R));report['resultant_degree']=int(s.degree(R,U));report['bidegrees']=[list(map(int,s.Poly(p,U,V).degree_list())) for p in (f,g)]
        else:
            H=s.Poly(f,U).gcd(s.Poly(g,U));report['univariate_gcd']=str(H.as_expr())
            if H.is_zero:raise NotImplementedError('An endpoint stationary component needs separate treatment')
            factors=s.factor_list(H.as_expr())[1]
        reports=[]
        for h,multiplicity in factors:
            hp=s.Poly(h,U);count=int(hp.count_roots(-s.oo,s.oo));row={'h':str(h),'multiplicity':int(multiplicity),'real_projection_roots':count}
            if not count:reports.append(row);continue
            F=Fiber(h)
            if chart=='finite':G=F.polynomial(f).gcd(F.polynomial(g)).sqf_part().monic()
            else:G=s.Poly.from_dict({(1,):F.K.one},V,domain=F.K) # artificial v=0; U is endpoint coordinate
            forbidden=[('normalization',s.fraction(geo['norm'])[0]),('phase_chart',1+V**2)]
            forbidden += [(f'current{i+1}',s.fraction(x)[0]) for i,x in enumerate(geo['h'])]
            if sector=='charged':forbidden += [('kernel_normalization',s.fraction(geo['W'])[0])]
            if chart=='finite':G,removed=remove_zero_factors(G,[(name,F.polynomial(pol)) for name,pol in forbidden])
            else:
                removed=[]
                if any(not F.evaluate_poly(pol,F.K.zero) for _,pol in forbidden):G=s.Poly(1,V,domain=F.K);removed=[('normalization_or_current',1)]
            row['removed']=removed;row['surviving_fiber']=str(F.poly_expr(G));reports.append(row)
            if G.degree():out['candidates']+=evaluate_candidates(h,F,G,geo,sector,chart)
        report['fibers']=reports;out['charts'][chart]=report
    return out

def algebraic_sign(x):
    x=s.simplify(x)
    if x==0:return 0
    if x.is_positive:return 1
    if x.is_negative:return -1
    y=s.sign(x)
    if y in (-1,0,1):return int(y)
    raise ArithmeticError('Could not certify the sign of an algebraic sparse value')

def sparse(data):
    a,b,c,dr,di=[[qinput(v) for v in data[k]]for k in ('a','b','c','d_real','d_imag')];A=s.diag(*a);out=[]
    for e,(i,j) in enumerate(E):A[i,j]=A[j,i]=b[e]
    def accept(label,B,support):
        n=len(support);one=s.ones(n,1);T=B.row_join(-one).col_join(one.T.row_join(s.zeros(1)));det=s.factor(T.det())
        if det==0:out.append({'label':label,'singular':'covered by boundary descent; not asserted absent'});return
        sol=T.inv()*s.Matrix([0]*n+[1]);r=[s.simplify(x)for x in sol[:-1]];mu=s.simplify(sol[-1]);physical=all(algebraic_sign(x)>0 for x in r);row={'label':label,'support':list(support),'r':list(map(str,r)),'value':str(s.simplify(mu/2)),'physical':physical,'value_sign':algebraic_sign(mu),'display_value':float(mu/2)};out.append(row)
    for i in range(3):accept(f'singleton_{i+1}',s.Matrix([[a[i]]]),(i,))
    for e,(i,j) in enumerate(E):
        mag=s.sqrt(dr[e]**2+di[e]**2)
        for sign in (-1,1):
            B=A.extract([i,j],[i,j]);B[0,1]=B[1,0]=b[e]+c[e]+sign*mag;accept(f'neutral_pair_{i+1}{j+1}_{sign}',B,(i,j))
        accept(f'orthogonal_pair_{i+1}{j+1}',A.extract([i,j],[i,j]),(i,j))
        for sign in (-1,1):
            B=A.copy();B[i,j]=B[j,i]=b[e]+c[e]+sign*mag;accept(f'split_{i+1}{j+1}_{sign}',B,(0,1,2))
    return out

def classify(data):
    out={'couplings':data,'scope':'Exact on the successfully audited finite nonreal-cycle, nonsingular-angular-block input; not an all-coupling classifier.','sparse':sparse(data),'sectors':{}}
    for sector in ('neutral','charged'):out['sectors'][sector]=solve_sector(data,sector)
    signs=[r['value_sign']for r in out['sparse']if r.get('physical')]+[r['signs']['mu']for sec in out['sectors'].values()for r in sec['candidates']if r['physical']]
    out['quartic_status']='unstable' if min(signs)<0 else ('nonnegative_with_flat_direction' if min(signs)==0 else 'strictly_positive');return out
if __name__=='__main__':
    path=Path(sys.argv[1]);data=json.loads(path.read_text());out=classify(data);target=Path(sys.argv[2])if len(sys.argv)>2 else path.with_suffix('.certificate.json');target.write_text(json.dumps(out,indent=2)+'\n');print(out['quartic_status']);print({k:(len(v['candidates']),sum(r['physical']for r in v['candidates']))for k,v in out['sectors'].items()})
