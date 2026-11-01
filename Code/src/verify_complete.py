"""Fail-closed audit of a decision, including finite-sample root completeness.

No classifier, geometry-constructor, or elimination-solver imports are used.
SymPy remains the shared arithmetic backend; this is not formal verification.
"""
from __future__ import annotations
import json
from pathlib import Path
import sympy as s
from sympy.polys.domains import QQ
from verify_sample import verify_sample
from verify_barrier import verify as verify_barrier
U,V=s.symbols('u v'); EDGES=((0,1),(1,2),(2,0)); OPP=(1,2,0)
KEYS=('a','b','c','d_real','d_imag'); ANCHOR=((5,6,7),(0,0,0),(2,3,4),(1,2,-1),(1,-1,2))

def canonical(f):
    p=s.Poly(s.cancel(f),U,V,domain=QQ)
    return p.monic().as_expr() if not p.is_zero else s.S.Zero

def direct_geometry(data,sector):
    a,b,c,dr,di=[[s.Rational(z) for z in data[k]] for k in KEYS]
    co=((1-U**2)/(1+U**2),((1-U**2)*(1-V**2)+4*U*V)/((1+U**2)*(1+V**2)),(1-V**2)/(1+V**2))
    si=(-2*U/(1+U**2),(2*U*(1-V**2)-2*V*(1-U**2))/((1+U**2)*(1+V**2)),2*V/(1+V**2))
    q=[dr[e]*co[e]-di[e]*si[e] for e in range(3)];w=[dr[e]*si[e]+di[e]*co[e] for e in range(3)];A=s.diag(*a)
    for e,(i,j) in enumerate(EDGES): A[i,j]=A[j,i]=b[e]
    if sector=='neutral':
        h=s.Matrix([w[e] for e in OPP]); B=A.copy()
        for e,(i,j) in enumerate(EDGES): B[i,j]=B[j,i]=b[e]+c[e]+q[e]
        m=B*h; norm=sum(h); W=norm
    elif sector=='charged':
        delta=[c[e]**2-dr[e]**2-di[e]**2 for e in range(3)]
        if any(z==0 for z in delta): raise AssertionError('Singular sample angular block')
        h=s.Matrix([-w[e]/delta[e] for e in OPP]); M=s.zeros(3)
        for e,(i,j) in enumerate(EDGES): M[i,j]=M[j,i]=(c[e]-q[e])/delta[e]
        k=(M*h).applyfunc(s.cancel); m=A*k+h; norm=sum(k); W=sum(h)
    else: raise AssertionError('Unknown sector')
    D=(1+U**2)*(1+V**2);f,g=[canonical(s.cancel(D*(m[i]-m[2]))) for i in (0,1)]
    bad=[('u_chart',1+U**2),('v_chart',1+V**2),('normalization',s.fraction(s.cancel(norm))[0])]
    bad += [(f'current{i}',s.fraction(s.cancel(z))[0]) for i,z in enumerate(h)]
    if sector=='charged': bad.append(('kernel_norm',s.fraction(s.cancel(W))[0]))
    return f,g,bad

def audit_endpoints(f,g):
    for x,y in ((U,V),(V,U)):
        F=s.Poly(s.Poly(f,x).nth(2),y,domain=QQ);G=s.Poly(s.Poly(g,x).nth(2),y,domain=QQ);H=F.gcd(G)
        if H.is_zero or (H.degree()>0 and H.count_roots(-s.oo,s.oo)>0): raise AssertionError('Unaccounted real endpoint')
    if s.Poly(f,U,V).coeff_monomial(U**2*V**2)==0 and s.Poly(g,U,V).coeff_monomial(U**2*V**2)==0: raise AssertionError('Unaccounted endpoint corner')

def audit_sector(data,sector,record):
    f,g,bad=direct_geometry(data,sector); report=record['charts']['finite']; audit_endpoints(f,g)
    if s.Poly(s.gcd(f,g),U,V).total_degree()>0: raise AssertionError('Nonfinite phase system')
    if [canonical(s.sympify(z,locals={'u':U,'v':V})) for z in report['equations']] != [f,g]: raise AssertionError('Wrong reported phase equations')
    k=report['projection_shear']
    if isinstance(k,bool) or not isinstance(k,int) or k<0: raise AssertionError('Invalid shear')
    F=s.expand(f.subs(U,U-k*V));G=s.expand(g.subs(U,U-k*V));R=s.Poly(s.resultant(F,G,V),U,domain=QQ)
    if R.is_zero: raise AssertionError('Zero eliminant')
    stored=s.Poly(s.sympify(report['resultant'],locals={'u':U}),U,domain=QQ)
    if stored.monic()!=R.monic(): raise AssertionError('Wrong eliminant')
    records=record['candidates']
    if any(r['chart']!='finite' or r['projection_shear']!=k for r in records): raise AssertionError('Unexpected chart or shear')
    seen=set(); expected_roots=0; factors=0; excluded=0
    for h,mul in s.factor_list(R.as_expr(),U)[1]:
        hp=s.Poly(h,U,domain=QQ).monic();nroots=int(hp.count_roots(-s.oo,s.oo));factors+=1
        matching=[(j,r) for j,r in enumerate(records) if s.Poly(s.sympify(r['h'],locals={'u':U}),U,domain=QQ).monic()==hp]
        if not nroots:
            if matching: raise AssertionError('Record for a nonreal projection factor')
            continue
        K=QQ.alg_field_from_poly(hp,alias='audit_alpha')
        def cv(poly):
            pp=s.Poly(poly,U,domain=QQ).rem(hp)
            return K.new([QQ.from_sympy(c) for c in pp.all_coeffs()])
        def pv(poly): return s.Poly.from_dict({powers:cv(c) for powers,c in s.Poly(poly,V).terms()},V,domain=K)
        fiber=pv(F).gcd(pv(G)).sqf_part().monic()
        for name,pol in bad:
            bp=pv(s.expand(pol.subs(U,U-k*V)))
            while fiber.degree()>0:
                common=fiber.gcd(bp)
                if common.degree()==0: break
                fiber=fiber.exquo(common)
        if fiber.degree()==0:
            if matching: raise AssertionError('Forbidden fiber has candidate records')
            excluded+=1;continue
        if fiber.degree()!=1: raise AssertionError('Unresolved nonlinear fiber')
        if len(matching)!=nroots: raise AssertionError('Missing or duplicated algebraic candidate')
        wanted=-fiber.rep.to_list()[1]/fiber.rep.to_list()[0]; boxes=[]
        for j,row in matching:
            if j in seen: raise AssertionError('Candidate reused')
            seen.add(j);lo,hi=[s.Rational(x) for x in row['u_interval']]
            if lo>hi or (lo==hi and hp.eval(lo)!=0) or (lo<hi and hp.count_roots(lo,hi)!=1): raise AssertionError('Invalid root interval')
            boxes.append((lo,hi))
            if cv(s.sympify(row['v_polynomial'],locals={'u':U}))!=wanted: raise AssertionError('Incorrect lifted phase coordinate')
            vals=row['field_values']
            if cv(s.sympify(vals['original_v'],locals={'u':U}))!=wanted: raise AssertionError('Wrong original v')
            if cv(s.sympify(vals['original_u'],locals={'u':U}))!=cv(U)-k*wanted: raise AssertionError('Wrong original u')
        boxes.sort()
        if any(boxes[i][1]>=boxes[i+1][0] for i in range(len(boxes)-1)): raise AssertionError('Duplicate or overlapping projection roots')
        expected_roots+=nroots
    if len(seen)!=len(records): raise AssertionError('Extraneous unaccounted candidate')
    return {'verified':True,'projection_shear':k,'eliminant_degree':R.degree(),'factor_count':factors,'excluded_real_fibers':excluded,'real_candidates_exhausted':expected_roots,'endpoint_exclusion_verified':True}

def audit_sparse_coverage(sample):
    expected={f'singleton_{i+1}':(i,) for i in range(3)}
    for i,j in EDGES:
        expected[f'orthogonal_pair_{i+1}{j+1}']=(i,j)
        for sign in (-1,1):expected[f'neutral_pair_{i+1}{j+1}_{sign}']=(i,j);expected[f'split_{i+1}{j+1}_{sign}']=(0,1,2)
    rows=sample['sparse'];labels=[r['label'] for r in rows]
    if len(labels)!=18 or set(labels)!=set(expected):raise AssertionError('Missing or duplicated sparse branch')
    for row in rows:
        if 'singular' in row:raise AssertionError('Singular sparse system contradicts sample regularity')
        if tuple(row['support'])!=expected[row['label']]:raise AssertionError('Wrong sparse support')

def audit_sample(sample):
    audit_sparse_coverage(sample);sign_report=verify_sample(sample);out={'candidate_sign_audit':sign_report,'coverage':{}}
    if sample['complete_enumeration']:
        if set(sample['sectors'])!={'neutral','charged'}:raise AssertionError('Missing entire full-support sector')
        for sector in ('neutral','charged'):out['coverage'][sector]=audit_sector(sample['couplings'],sector,sample['sectors'][sector])
    else:
        witnesses=[r for r in sample['sparse'] if r.get('physical') and r['value_sign']<0]
        if sample['quartic_status']!='unstable' or not witnesses:raise AssertionError('Incomplete enumeration lacks valid negative witness')
        out['early_negative_witness']=witnesses[0]['label']
    return out

def audit_decision(record,check_barrier=True):
    out={}
    if check_barrier:out['barrier']=verify_barrier(record['barrier'])
    if record['couplings']!=record['barrier']['couplings']:raise AssertionError('Original/barrier input mismatch')
    eta=s.Rational(record['barrier']['sample_epsilon'])
    if eta<=0 or set(record['samples'])!={'1','-1'}:raise AssertionError('Invalid sample pair')
    out['samples']={}
    for side in (1,-1):
        sample=record['samples'][str(side)]
        for key,q in zip(KEYS,ANCHOR):
            expected=[s.Rational(x)+side*eta*y for x,y in zip(record['couplings'][key],q)]
            if expected!=[s.Rational(x) for x in sample['couplings'][key]]:raise AssertionError('Sample is not the declared deformation')
        out['samples'][str(side)]=audit_sample(sample)
    pair=tuple(record['samples'][str(side)]['quartic_status'] for side in (1,-1))
    table={('strictly_positive','strictly_positive'):'strictly_positive',('strictly_positive','unstable'):'nonnegative_with_flat_direction',('unstable','unstable'):'unstable'}
    if table.get(pair)!=record['quartic_status']:raise AssertionError('Invalid three-way decision')
    out.update(verified=True,status=record['quartic_status'],scope='Complete finite-sample enumeration and arithmetic audit, not independent specialist/formal verification.')
    return out

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('certificate');ap.add_argument('--output');args=ap.parse_args()
    out=audit_decision(json.loads(Path(args.certificate).read_text()));txt=json.dumps(out,indent=2)
    if args.output:Path(args.output).write_text(txt+'\n')
    print(txt)
