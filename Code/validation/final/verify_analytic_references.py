"""Exact reference identities independent of current equations and resultants."""
from pathlib import Path
import sys,json
import sympy as s
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(Path(__file__).parent))
from run_closeout_cases import cases
r1,r2,r3,x12,y12,x23,y23,x31,y31=s.symbols('r1 r2 r3 x12 y12 x23 y23 x31 y31', real=True)
r=s.Matrix([r1,r2,r3]); S=sum(r); xs=[x12,x23,x31];ys=[y12,y23,y31];edges=[(0,1),(1,2),(2,0)]
def poly(p):
    a,b,c,dr,di=[[s.Rational(x) for x in p[k]] for k in ('a','b','c','d_real','d_imag')]
    return sum(a[i]*r[i]**2/2 for i in range(3))+sum(b[e]*r[i]*r[j]+c[e]*(xs[e]**2+ys[e]**2)+dr[e]*(xs[e]**2-ys[e]**2)-2*di[e]*xs[e]*ys[e] for e,(i,j) in enumerate(edges))
f=cases();p=f['rotated_singular_flat']['input']
sos=2*(r1-r2/5)**2+3*(r2/5-r3)**2+5*(r3-r1)**2+(x12+2*y12)**2+2*(2*x23-3*y23)**2+3*(3*x31+y31)**2
assert s.expand(poly(p)-sos)==0
X=s.Matrix([[1,2-s.I,0],[0,0,1]]);K=s.conjugate(X.T)*X
subs={r1:K[0,0],r2:K[1,1],r3:K[2,2],x12:s.re(K[0,1]),y12:s.im(K[0,1]),x23:s.re(K[1,2]),y23:s.im(K[1,2]),x31:s.re(K[2,0]),y31:s.im(K[2,0])}
assert s.simplify(sos.subs(subs))==0
for name in ('rotated_singular_positive','rotated_singular_negative'):
    shift=s.Rational(f[name]['exact_minimum']);assert s.expand(poly(f[name]['input'])-sos-shift*S**2)==0
A=s.Matrix([[4,1,-1],[1,6,0],[-1,0,8]]); rr=s.Matrix([s.Rational(23,48),s.Rational(11,48),s.Rational(7,24)]); m=s.Rational(89,96)
assert sum(rr)==1 and A*rr==2*m*s.ones(3,1)
minors=[A[:i,:i].det() for i in (1,2,3)];assert all(v>0 for v in minors)
assert s.expand(poly(f['radial_nondiagonal_positive']['input'])-m*S**2-((r-S*rr).T*A*(r-S*rr))[0]/2)==0
for name in ('radial_nondiagonal_flat','radial_nondiagonal_negative'):
    target=s.Rational(f[name]['exact_minimum']);assert s.expand(poly(f[name]['input'])-poly(f['radial_nondiagonal_positive']['input'])-(target-m)*S**2)==0
out={'verified':True,'reference_cases':6,'method':'Exact expansion of independently chosen positive-square identities, an explicit physical zero, and Sylvester principal-minor test; no stationary solver or resultant imported.', 'rotated_singular':{'sum_of_squares_identity':True,'X_real':[[1,2,0],[0,0,1]],'X_imag':[[0,-1,0],[0,0,0]],'trace_gram':str(s.simplify(s.trace(K))),'rank_gram':K.rank(),'cycle_product':str(s.prod(s.Rational(a)+s.I*s.Rational(b) for a,b in zip(p['d_real'],p['d_imag'])) .expand()),'deltas':[str(s.Rational(c)**2-s.Rational(a)**2-s.Rational(b)**2) for c,a,b in zip(p['c'],p['d_real'],p['d_imag'])]},'radial':{'leading_principal_minors':[str(v) for v in minors],'minimizing_radii':[str(v) for v in rr],'normalized_minimum':str(m),'lower_bound_identity_verified':True},'reference_minima':{k:v['exact_minimum'] for k,v in f.items()}}
(ROOT/'validation/final/analytic_reference_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
