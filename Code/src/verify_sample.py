"""Audit recorded finite-sample candidates using exact quotient arithmetic.

Does not import the solvers. Checks isolating intervals, all recorded algebraic
signs, direct angular-matrix reconstruction, physicality, and sparse stationary
records. This checks returned candidates and interval coverage in the records;
it is NOT a separate proof of exhaustiveness of the polynomial elimination.
"""
from __future__ import annotations
import json
from pathlib import Path
import sympy as s
from sympy.polys.domains import QQ
u=s.Symbol('u');edges=((0,1),(1,2),(2,0));opp=(1,2,0)


def sign_algebraic(expr):
    expr=s.simplify(expr)
    if expr==0:return 0
    sign=s.sign(expr)
    if sign in (-1,1):return int(sign)
    raise ArithmeticError('Unresolved exact algebraic sign')


def interval_eval(p,lo,hi):
    a=b=s.S.Zero
    for c in p.all_coeffs():
        vals=(a*lo,a*hi,b*lo,b*hi);a=min(vals)+c;b=max(vals)+c
    return a,b


def verify_sample(sample):
    a,b,c,dr,di=[[s.Rational(x)for x in sample['couplings'][k]]for k in ('a','b','c','d_real','d_imag')];A=s.diag(*a);sparse_count=0;root_count=0;inequalities=0;physical=0
    for j,(i,k)in enumerate(edges):A[i,k]=A[k,i]=b[j]
    for row in sample['sparse']:
        if 'singular'in row:continue
        label=row['label'];support=row['support'];B=A.extract(support,support);r=s.Matrix([s.sympify(x)for x in row['r']]);value=s.sympify(row['value'])
        if label.startswith(('neutral_pair','split')):
            j=next(j for j,(i,k)in enumerate(edges)if str(i+1)+str(k+1)in label.rsplit('_',2)[1]);i,k=edges[j];sgn=int(label.rsplit('_',1)[1]);co=b[j]+c[j]+sgn*s.sqrt(dr[j]**2+di[j]**2)
            if label.startswith('neutral_pair'):B[0,1]=B[1,0]=co
            else:B[i,k]=B[k,i]=co
        if any(s.simplify(x)!=0 for x in B*r-2*value*s.ones(len(support),1)) or s.simplify(sum(r)-1)!=0:raise AssertionError('Sparse stationary identity failed')
        good=all(sign_algebraic(x)>0 for x in r)
        if good!=row['physical'] or sign_algebraic(value)!=row['value_sign']:raise AssertionError('Sparse sign mismatch')
        sparse_count+=1
    for sector,sec in sample['sectors'].items():
        for rec in sec['candidates']:
            if rec['chart']!='finite':raise NotImplementedError('This checkpoint barrier excludes endpoint roots')
            hp=s.Poly(s.sympify(rec['h'],locals={'u':u}),u,domain=QQ);lo,hi=map(s.Rational,rec['u_interval'])
            if lo==hi:
                if hp.eval(lo)!=0:raise AssertionError('Invalid exact root')
            elif hp.count_roots(lo,hi)!=1:raise AssertionError('Invalid isolating interval')
            K=QQ.alg_field_from_poly(hp,alias='a');x=K.new([QQ.one,QQ.zero]);zero=K.zero;one=K.one
            def cv(expr):
                p=s.Poly(s.sympify(expr,locals={'u':u}),u,domain=QQ).rem(hp);return K.new([QQ.from_sympy(z)for z in p.all_coeffs()])
            vals={key:cv(expr)for key,expr in rec['field_values'].items()}
            for key,expr in rec['field_values'].items():
                p=s.Poly(s.sympify(expr,locals={'u':u}),u,domain=QQ).rem(hp);want=rec['signs'].get(key)
                if want is None:continue
                if p.is_zero:actual=0
                else:
                    low,high=interval_eval(p,lo,hi)
                    if low>0:actual=1
                    elif high<0:actual=-1
                    else:raise AssertionError('Stored interval fails to certify '+key)
                if want!=actual:raise AssertionError('Incorrect algebraic sign')
                inequalities+=1
            vv=vals.get('original_v',vals.get('v',zero));uu=vals.get('original_u',x-rec.get('projection_shear',0)*vv)
            cs=lambda t:(one-t*t)/(one+t*t);sn=lambda t:2*t/(one+t*t)
            cx,sx,cy,sy=cs(uu),sn(uu),cs(vv),sn(vv);co=(cx,cx*cy+sx*sy,cy);si=(-sx,sx*cy-cx*sy,sy)
            q=[dr[j]*co[j]-di[j]*si[j]for j in range(3)];w=[dr[j]*si[j]+di[j]*co[j]for j in range(3)]
            if sector=='neutral':
                h=[w[j]for j in opp];W=sum(h,zero);rr=[x/W for x in h];m=[a[i]*rr[i]for i in range(3)]
                for j,(i,k)in enumerate(edges):m[i]+=(b[j]+c[j]+q[j])*rr[k];m[k]+=(b[j]+c[j]+q[j])*rr[i]
                expected={'r'+str(i+1):rr[i]for i in range(3)};expected['mu']=m[0]
            else:
                delta=[c[j]**2-dr[j]**2-di[j]**2 for j in range(3)];h=[-w[j]/delta[j]for j in opp];kk=[zero,zero,zero]
                for j,(i,k)in enumerate(edges):coef=(c[j]-q[j])/delta[j];kk[i]+=coef*h[k];kk[k]+=coef*h[i]
                C=sum(kk,zero);rr=[x/C for x in kk];zz=[-x/C for x in h];nu=sum(zz,zero);ss=[z/nu for z in zz]
                m=[sum((A[i,j]*rr[j]for j in range(3)),zero)-zz[i]for i in range(3)]
                margin=one-sum((x*x for x in rr),zero)
                for j,(i,k)in enumerate(edges):margin-=2*zz[i]*zz[k]*(c[j]**2+dr[j]**2+di[j]**2-2*c[j]*q[j])/delta[j]**2
                expected={'r'+str(i+1):rr[i]for i in range(3)};expected.update({'s'+str(i+1):ss[i]for i in range(3)});expected.update(mu=m[0],nu=nu,psd_margin=margin)
            if any(z!=m[0]for z in m)or sum(rr,zero)!=one:raise AssertionError('Direct stationary matrix reconstruction failed')
            if any(vals[key]!=value for key,value in expected.items()):raise AssertionError('Recorded algebraic quantity mismatch')
            good=all(rec['signs']['r'+str(i)]>0 for i in (1,2,3))
            if sector=='charged':good=good and all(rec['signs']['s'+str(i)]>0 for i in (1,2,3))and rec['signs']['psd_margin']>=0
            if good!=rec['physical']:raise AssertionError('Physicality mismatch')
            root_count+=1;physical+=int(good)
    signs=[r['value_sign']for r in sample['sparse']if r.get('physical')]+[r['signs']['mu']for sec in sample['sectors'].values()for r in sec['candidates']if r['physical']]
    inferred='unstable'if min(signs)<0 else ('strictly_positive'if min(signs)>0 else 'nonnegative_with_flat_direction')
    if inferred!=sample['quartic_status']:raise AssertionError('Sample status mismatch')
    if inferred!='unstable'and not sample['complete_enumeration']:raise AssertionError('Positive sample needs complete enumeration')
    return {'verified':True,'sparse_records':sparse_count,'algebraic_roots_checked':root_count,'physical_full_support_roots':physical,'exact_algebraic_signs_checked':inequalities,'status':inferred,'scope':'Candidate and exact-sign audit, not a new elimination-completeness proof.'}


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('certificate');ap.add_argument('--output');args=ap.parse_args();raw=json.loads(Path(args.certificate).read_text());out={side:verify_sample(v)for side,v in raw['samples'].items()};txt=json.dumps(out,indent=2);print(txt)
    if args.output:Path(args.output).write_text(txt+'\n')
