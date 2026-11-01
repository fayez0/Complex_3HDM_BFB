"""Checkpoint-03 finite-current solver with exact separating projections.

Original checkpoint-02 code is imported unchanged. Nonlinear lifting fibers are
resolved by trying t=u+k*v, with exact substitution and sign certification.
Continuous components and singular angular blocks still raise explicit errors.
"""
from __future__ import annotations
import json,sys
from pathlib import Path
import sympy as s
from exact_current import geometry,evaluate_candidates,sparse
from exact_phase_solver import U,V,Fiber,remove_zero_factors


def _shear_geometry(geo,k):
    out=dict(geo)
    for name in ('f','g','mu','norm','W'):
        out[name]=s.expand(geo[name].subs(U,U-k*V)) if name in ('f','g') else geo[name].subs(U,U-k*V)
    for name in ('r','z','h','q'): out[name]=[z.subs(U,U-k*V) for z in geo[name]]
    return out


def solve_finite(geo,sector,max_shears=121):
    if s.Poly(s.gcd(geo['f'],geo['g']),U,V).total_degree()>0: raise NotImplementedError('Continuous phase component: use the universal deformation theorem, not deletion')
    attempts=[]
    for k in range(max_shears):
        ge=_shear_geometry(geo,k);f,g=ge['f'],ge['g'];R=s.resultant(f,g,V)
        if R==0: attempts.append({'k':k,'status':'zero resultant'});continue
        candidates=[];fibers=[];nonlinear=False
        for h,multiplicity in s.factor_list(R)[1]:
            hp=s.Poly(h,U);count=int(hp.count_roots(-s.oo,s.oo));row={'h':str(h),'multiplicity':int(multiplicity),'real_projection_roots':count}
            if not count: fibers.append(row);continue
            F=Fiber(h);G=F.polynomial(f).gcd(F.polynomial(g)).sqf_part().monic()
            bad=[('original_u_chart',1+(U-k*V)**2),('original_v_chart',1+V**2),('normalization',s.fraction(s.cancel(ge['norm']))[0])]
            bad += [(f'current{i+1}',s.fraction(s.cancel(x))[0]) for i,x in enumerate(ge['h'])]
            if sector=='charged':bad += [('kernel_normalization',s.fraction(s.cancel(ge['W']))[0])]
            G,removed=remove_zero_factors(G,[(name,F.polynomial(pol)) for name,pol in bad]);row.update(removed=removed,surviving_fiber=str(F.poly_expr(G)));fibers.append(row)
            if G.degree()>1: nonlinear=True;row['reason']='nonseparating projection';break
            if G.degree()==1:
                rows=evaluate_candidates(h,F,G,ge,sector,'finite');vv=-G.rep.to_list()[1]/G.rep.to_list()[0];uu=F.a-k*vv
                for record in rows:
                    box=tuple(s.Rational(x) for x in record['u_interval']);_,iv,_=F.signs({'original_u':uu,'original_v':vv},box)
                    record['projection_shear']=k;record['projection_definition']=f't = original_u + {k}*original_v'
                    record['field_values']['original_u']=str(F.as_expr(uu));record['field_values']['original_v']=str(F.as_expr(vv))
                    record['display_tangents']=[float(sum(iv['original_u'])/2),float(sum(iv['original_v'])/2)]
                candidates+=rows
        if nonlinear: attempts.append({'k':k,'status':'nonlinear surviving fiber','fiber_audit':fibers});continue
        attempts.append({'k':k,'status':'all real projection fibers exhausted'})
        return {'equations':[str(geo['f']),str(geo['g'])],'projection_shear':k,'resultant':str(s.factor(R)),'resultant_degree':int(s.degree(R,U)),'fibers':fibers,'projection_attempts':attempts},candidates
    raise ArithmeticError('Separating-projection search exceeded its declared resource budget')


def solve_sector(data,sector,max_shears=121):
    out={'sector':sector,'charts':{},'candidates':[]};geo=geometry(data,sector,'finite');report,rows=solve_finite(geo,sector,max_shears);out['charts']['finite']=report;out['candidates']+=rows
    for chart in ('u_infinity','v_infinity','corner'):
        geo=geometry(data,sector,chart);f,g=geo['f'],geo['g'];report={'equations':[str(f),str(g)]}
        if chart=='corner':
            report['stationary']=f==g==0
            if report['stationary']:
                if geo['norm']==0 or geo['W']==0 or any(x==0 for x in geo['h']):report['excluded']='zero normalization or current'
                else:
                    F=Fiber(U);G=s.Poly.from_dict({(1,):F.K.one},V,domain=F.K);out['candidates']+=evaluate_candidates(U,F,G,geo,sector,chart)
            out['charts'][chart]=report;continue
        H=s.Poly(f,U).gcd(s.Poly(g,U));report['univariate_gcd']=str(H.as_expr())
        if H.is_zero:raise NotImplementedError('An endpoint stationary component requires deformation or direct component treatment')
        fibers=[]
        for h,multiplicity in s.factor_list(H.as_expr())[1]:
            count=int(s.Poly(h,U).count_roots(-s.oo,s.oo));row={'h':str(h),'multiplicity':int(multiplicity),'real_projection_roots':count}
            if not count:fibers.append(row);continue
            F=Fiber(h);bad=[s.fraction(geo['norm'])[0],1+U**2]+[s.fraction(x)[0] for x in geo['h']]
            if sector=='charged':bad+=[s.fraction(geo['W'])[0]]
            if any(not F.evaluate_poly(z,F.K.zero) for z in bad):row['excluded']='zero normalization or current'
            else:
                G=s.Poly.from_dict({(1,):F.K.one},V,domain=F.K);out['candidates']+=evaluate_candidates(h,F,G,geo,sector,chart)
            fibers.append(row)
        report['fibers']=fibers;out['charts'][chart]=report
    return out


def classify(data,max_shears=121):
    out={'couplings':data,'scope':'Exact finite nonreal-cycle classifier with nonsingular angular blocks; separating projections supported; arbitrary infinitesimal classification not implemented.','sparse':sparse(data),'sectors':{}}
    for sector in ('neutral','charged'):out['sectors'][sector]=solve_sector(data,sector,max_shears)
    signs=[r['value_sign'] for r in out['sparse'] if r.get('physical')]+[r['signs']['mu'] for sec in out['sectors'].values() for r in sec['candidates'] if r['physical']]
    out['quartic_status']='unstable' if min(signs)<0 else ('nonnegative_with_flat_direction' if min(signs)==0 else 'strictly_positive');return out

if __name__=='__main__':
    data=json.loads(Path(sys.argv[1]).read_text());out=classify(data);target=Path(sys.argv[2]);target.write_text(json.dumps(out,indent=2)+'\n');print(out['quartic_status'])
