"""Fast API and exact-algebra regressions. Full classifications live in validation/final."""
from pathlib import Path
from fractions import Fraction
import copy, importlib.util, json, subprocess, sys, tempfile, unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'src'))
from weinberg_bfb import validate_coefficients,load_json,audit_any
BASE={'a':[5,6,7],'b':[0,0,0],'c':[2,3,4],'d_real':[1,2,-1],'d_imag':[1,-1,2]}

class InputContract(unittest.TestCase):
    def test_rational(self):
        p=copy.deepcopy(BASE);p['a']=['1/3','0.125','1e-30'];self.assertEqual(validate_coefficients(p)['a'],['1/3','1/8','1/'+str(10**30)])
    def test_reject_float(self):
        p=copy.deepcopy(BASE);p['a'][0]=0.5
        with self.assertRaises(ValueError):validate_coefficients(p)
    def test_reject_bool(self):
        p=copy.deepcopy(BASE);p['a'][0]=True
        with self.assertRaises(ValueError):validate_coefficients(p)
    def test_missing_and_extra(self):
        for p in ({k:v for k,v in BASE.items() if k!='c'},dict(BASE,extra=[])):
            with self.assertRaises(ValueError):validate_coefficients(p)
    def test_wrong_shapes(self):
        for row in ([1,2], [1,2,3,4], '123', (1,2,3)):
            p=copy.deepcopy(BASE);p['a']=row
            with self.assertRaises(ValueError):validate_coefficients(p)
    def test_nonfinite_and_expressions(self):
        for bad in ('NaN','Infinity','1/0','sqrt(2)','__import__("os").getcwd()'):
            p=copy.deepcopy(BASE);p['a'][0]=bad
            with self.assertRaises(ValueError):validate_coefficients(p)
    def test_duplicate_json_keys(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/'bad.json';f.write_text('{"a":1,"a":2}')
            with self.assertRaises(ValueError):load_json(f)

class ExactAlgebra(unittest.TestCase):
    def test_nonradical_fiber(self):
        import sympy as s
        from exact_phase_solver import U,V,Fiber
        F=Fiber(U*U-2);h=F.polynomial((V-U)**2).sqf_part().monic()
        self.assertEqual(h.degree(),1);value=-h.rep.to_list()[1]/h.rep.to_list()[0]
        self.assertEqual(value,F.a);self.assertFalse(F.evaluate_poly((V-U)**2,value))
    def test_forbidden_fiber_saturation(self):
        from exact_phase_solver import U,V,Fiber,remove_zero_factors
        F=Fiber(U-1);h=F.polynomial((V-1)**3*(V+2)**2).sqf_part()
        reduced,_=remove_zero_factors(h,[('forbidden',F.polynomial(V-1))])
        self.assertEqual(reduced.degree(),1);self.assertEqual(-reduced.rep.to_list()[1]/reduced.rep.to_list()[0],F.K.convert(-2))
    def test_dominance_fallback(self):
        import sympy as s
        from robust_sign import dominance_radius
        from zero_barrier import T,no_root_test
        p=s.Poly(T+s.Rational(1,10**100),T);eta,record=dominance_radius({'small':p})
        self.assertGreater(eta,0);self.assertLess(eta,s.Rational(1,10**100))
        for side in (-1,1):self.assertTrue(no_root_test(p,eta,side)[0])
    def test_anchor_degrees(self):
        from zero_barrier import anchor_zero_certificate
        rec=anchor_zero_certificate();self.assertEqual(rec['factor_count'],35);self.assertEqual(rec['product_degree'],517)
        self.assertTrue(all(Fraction(r['value'])!=0 for r in rec['factors'].values()))

class InterfaceExecution(unittest.TestCase):
    def run_cli(self,command,*args,stdlib=False):
        return subprocess.run([sys.executable]+(['-S'] if stdlib else [])+[str(ROOT/'weinberg_bfb.py'),command,*map(str,args)],text=True,capture_output=True,timeout=90)
    def test_stdlib_witness(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'out.json';p=self.run_cli('verify-witness',ROOT/'examples/Bminus_small_witness.json','--output',out,stdlib=True)
            self.assertEqual(p.returncode,0,p.stderr);self.assertTrue(load_json(out)['verified'])
    def test_existing_destination_protected(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'out.json';out.write_text('unchanged');p=self.run_cli('verify-witness',ROOT/'examples/Bminus_small_witness.json','--output',out)
            self.assertNotEqual(p.returncode,0);self.assertEqual(out.read_text(),'unchanged')
    def test_failure_has_no_sign(self):
        with tempfile.TemporaryDirectory() as d:
            inp=Path(d)/'in.json';out=Path(d)/'out.json';p=copy.deepcopy(BASE);p['a'][0]=0.5;inp.write_text(json.dumps(p))
            run=self.run_cli('classify',inp,'--output',out);self.assertNotEqual(run.returncode,0);record=load_json(out)
            self.assertFalse(record['completed']);self.assertNotIn('quartic_status',record)
    def test_stale_output_replaced_on_forced_failure(self):
        with tempfile.TemporaryDirectory() as d:
            inp=Path(d)/'in.json';out=Path(d)/'out.json';inp.write_text('{}');out.write_text('{"quartic_status":"strictly_positive"}')
            run=self.run_cli('classify',inp,'--output',out,'--force');self.assertNotEqual(run.returncode,0);self.assertNotIn('quartic_status',load_json(out))
    def test_optimization_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'out.json';run=subprocess.run([sys.executable,'-O',str(ROOT/'weinberg_bfb.py'),'verify-witness',str(ROOT/'examples/Bminus_small_witness.json'),'--output',str(out)],capture_output=True,text=True)
            self.assertNotEqual(run.returncode,0);self.assertFalse(load_json(out)['completed'])
    def test_complete_stored_negative_audit(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'audit.json';run=self.run_cli('verify',ROOT/'validation/extreme_scale/negative_normalized_classification.json','--output',out)
            self.assertEqual(run.returncode,0,run.stderr);self.assertTrue(load_json(out)['verified'])

if __name__=='__main__':unittest.main()
