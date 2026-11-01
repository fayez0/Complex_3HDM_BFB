"""New exact classifications with analytic reference answers, not historical reruns."""
from pathlib import Path
import json, sys, time, subprocess, traceback
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'validation/final'; sys.path.insert(0,str(ROOT/'src'))
KEYS=('a','b','c','d_real','d_imag')

def cases():
    from fractions import Fraction as Q
    sos={'a':[14,'2/5',16],'b':['-4/5','-6/5',-10],'c':['5/2',13,15],'d_real':['-3/2',-5,12],'d_imag':[-2,12,-9]}
    def shift(p,k):return {key:[str(Q(x)+(2*Q(k) if key in ('a','b') else 0)) for x in p[key]] for key in KEYS}
    radial={'a':[4,6,8],'b':[1,0,-1],'c':[0,0,0],'d_real':[0,0,0],'d_imag':[0,0,0]}
    return {
      'rotated_singular_flat':{'input':sos,'expected':'nonnegative_with_flat_direction','exact_minimum':'0','reference':'Three positive radial squares plus three rotated real-overlap squares; rational rank-two zero fields.'},
      'rotated_singular_positive':{'input':shift(sos,'1/10000'),'expected':'strictly_positive','exact_minimum':'1/10000','reference':'Exact constant shift of a manifest sum of squares with an attaining zero.'},
      'rotated_singular_negative':{'input':shift(sos,'-1/10000'),'expected':'unstable','exact_minimum':'-1/10000','reference':'Exact constant shift; the same rational fields attain the exact negative minimum.'},
      'radial_nondiagonal_positive':{'input':radial,'expected':'strictly_positive','exact_minimum':'89/96','reference':'Positive-definite radial matrix, positive stationary radii (23/48,11/48,7/24). Angular variables are unrestricted.'},
      'radial_nondiagonal_flat':{'input':shift(radial,'-89/96'),'expected':'nonnegative_with_flat_direction','exact_minimum':'0','reference':'Subtract the exact radial minimum times trace squared; a continuum of angular minima remains.'},
      'radial_nondiagonal_negative':{'input':shift(radial,'-179/192'),'expected':'unstable','exact_minimum':'-1/192','reference':'Constant shift below the exact radial minimum.'}}

def worker(name):
    from scale_normalized_sign import classify,audit
    fixture=cases()[name]; start=time.monotonic(); folder=OUT/'cases';folder.mkdir(exist_ok=True)
    (folder/(name+'_input.json')).write_text(json.dumps(fixture['input'],indent=2)+'\n')
    try:
        d=classify(fixture['input']); (folder/(name+'_decision.json')).write_text(json.dumps(d,indent=2)+'\n');tc=time.monotonic()-start
        if d['quartic_status']!=fixture['expected']:raise AssertionError('Analytic reference mismatch')
        a=audit(d); (folder/(name+'_audit.json')).write_text(json.dumps(a,indent=2)+'\n')
        result={**fixture,'passed':bool(a['verified']),'status':d['quartic_status'],'classification_seconds':tc,'total_seconds':time.monotonic()-start,'original_eta':d['original_deformation_radius'],'complete_audit':a['verified']}
    except Exception as exc:result={**fixture,'passed':False,'error':repr(exc),'traceback':traceback.format_exc(),'total_seconds':time.monotonic()-start}
    (folder/(name+'_result.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    return 0 if result['passed'] else 1

def main():
    results={};(OUT/'analytic_inputs.json').write_text(json.dumps(cases(),indent=2)+'\n')
    for name,fixture in cases().items():
        print('START',name,flush=True);start=time.monotonic()
        try:
            done=subprocess.run([sys.executable,__file__,'--worker',name],capture_output=True,text=True,timeout=300)
            (OUT/(name+'.log')).write_text(done.stdout+done.stderr)
            target=OUT/'cases'/(name+'_result.json')
            results[name]=json.loads(target.read_text()) if target.exists() else {**fixture,'passed':False,'returncode':done.returncode,'error':'No complete result emitted'}
        except subprocess.TimeoutExpired:
            results[name]={**fixture,'passed':False,'error':'Resource timeout (300 s); no verdict claimed','total_seconds':time.monotonic()-start}
        (OUT/'analytic_results.json').write_text(json.dumps(results,indent=2)+'\n');print('DONE',name,results[name].get('passed'),round(time.monotonic()-start,3),flush=True)
    return 0 if all(r['passed'] for r in results.values()) else 1
if __name__=='__main__':sys.exit(worker(sys.argv[2]) if len(sys.argv)>1 and sys.argv[1]=='--worker' else main())
