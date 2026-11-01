"""Add original-potential rational witnesses or explicit positive margins.

All inherited classifiers are unchanged. Numerical reconstruction is used ONLY
to propose rational fields; exact standard-library arithmetic accepts a negative
witness. A positive margin requires the complete existing decision audit.
"""
from __future__ import annotations
from pathlib import Path
import argparse, json, math
import mpmath as mp
import sympy as s
from verify_rational_witness import evaluate,verify,KEYS,EDGES
U=s.Symbol('u')

def numeric(expr,dps,alpha=None):
    value=s.sympify(expr,locals={'u':U})
    if alpha is not None: value=value.subs(U,alpha)
    return mp.mpf(str(s.N(value,dps)))

def full_fields(row,data,bits):
    dps=max(60,int(bits*.302)+40);mp.mp.dps=dps
    p=s.Poly(s.sympify(row['h'],locals={'u':U}),U);box=tuple(s.Rational(z) for z in row['u_interval'])
    if box[0]!=box[1]: box=p.refine_root(*box,eps=s.Rational(1,10**(dps+20)))
    alpha=(box[0]+box[1])/2;vals={k:numeric(v,dps,alpha) for k,v in row['field_values'].items()}
    u=vals.get('original_u',numeric(str(alpha),dps));v=vals.get('original_v',vals.get('v',mp.mpf(0)))
    half=[mp.mpc(1),(1+1j*u)/mp.sqrt(1+u*u),(1+1j*v)/mp.sqrt(1+v*v)]
    r=[vals['r'+str(i)] for i in (1,2,3)]
    if row['sector']=='neutral': return [[mp.sqrt(r[i])*mp.conj(half[i]) for i in range(3)],[mp.mpc(0)]*3]
    nn=[mp.sqrt(vals['s'+str(i+1)])*half[i] for i in range(3)];nu=vals['nu'];gram=mp.matrix(3)
    for i in range(3): gram[i,i]=r[i]
    for e,(i,j) in enumerate(EDGES):
        c=numeric(data['c'][e],dps);d=mp.mpc(numeric(data['d_real'][e],dps),numeric(data['d_imag'][e],dps));delta=c*c-abs(d)**2
        gram[i,j]=nu*(c*nn[i]*mp.conj(nn[j])-mp.conj(d)*mp.conj(nn[i])*nn[j])/delta;gram[j,i]=mp.conj(gram[i,j])
    # Two-row Gram factorization. Rounding need not preserve the Gram exactly:
    # the final rational fields are physical by construction and checked anew.
    pivot=max(range(3),key=lambda i:r[i]);first=[gram[pivot,j]/mp.sqrt(r[pivot]) for j in range(3)]
    residual=mp.matrix([[gram[i,j]-mp.conj(first[i])*first[j] for j in range(3)] for i in range(3)])
    k=max(range(3),key=lambda i:mp.re(residual[i,i]));diag=mp.re(residual[k,k])
    second=[residual[k,j]/mp.sqrt(diag) for j in range(3)] if diag>mp.mpf(10)**(-dps//2) else [mp.mpc(0)]*3
    return [first,second]

def sparse_fields(row,data,bits):
    dps=max(60,int(bits*.302)+40);mp.mp.dps=dps;X=[[mp.mpc(0) for _ in range(3)] for _ in range(2)]
    label=row['label'];r=[numeric(z,dps) for z in row['r']];support=row['support']
    if label.startswith('singleton'): X[0][support[0]]=mp.sqrt(r[0]);return X
    if label.startswith('orthogonal'): X[0][support[0]]=mp.sqrt(r[0]);X[1][support[1]]=mp.sqrt(r[1]);return X
    pair=label.split('_')[-2];e=next(e for e,(i,j) in enumerate(EDGES) if f'{i+1}{j+1}'==pair);i,j=EDGES[e];sign=int(label.split('_')[-1])
    radii={idx:r[k] for k,idx in enumerate(support)};d=mp.mpc(numeric(data['d_real'][e],dps),numeric(data['d_imag'][e],dps));phase=mp.sqrt(sign*mp.conj(d)/abs(d))
    X[0][i]=mp.sqrt(radii[i]);X[0][j]=mp.sqrt(radii[j])*phase
    if label.startswith('split'):
        k=next(k for k in range(3) if k not in (i,j));X[1][k]=mp.sqrt(radii[k])
    return X

def rational_witness(decision,max_bits=512):
    if decision['quartic_status']!='unstable': raise ValueError('Need an unstable original input')
    sample=decision['samples']['1'];p=decision['couplings']
    candidates=[('sparse',row) for row in sample['sparse'] if row.get('physical') and row['value_sign']<0]
    candidates += [('full',row) for sec in sample['sectors'].values() for row in sec['candidates'] if row['physical'] and row['signs']['mu']<0]
    if not candidates: raise ValueError('The positive deformation has no negative physical candidate')
    candidates.sort(key=lambda item:item[1].get('display_value',0));attempts=[]
    for bits in (16,32,64,128,256,512):
        if bits>max_bits: break
        for kind,row in candidates:
            X=sparse_fields(row,sample['couplings'],bits) if kind=='sparse' else full_fields(row,sample['couplings'],bits)
            den=2**bits;R=[[str(s.Rational(int(mp.nint(mp.re(x)*den)),den)) for x in line] for line in X];I=[[str(s.Rational(int(mp.nint(mp.im(x)*den)),den)) for x in line] for line in X]
            value,norm,_=evaluate(p,R,I);attempts.append({'bits':bits,'kind':kind,'negative':value<0})
            if value<0:
                out={'couplings':p,'fields_real':R,'fields_imag':I,'quartic_value':str(value),'field_norm':str(norm),'normalized_value':str(value/(norm*norm)),'componentwise_robustness_radius':str(-3*value/(4*norm*norm)),'source':'Physical candidate of positive deformation, rounded to rational ORIGINAL Higgs fields','rounding_bits':bits,'attempts':attempts}
                out['standalone_verification']=verify(out);return out
    raise ArithmeticError('Could not produce a rational witness within the stated precision budget; no physical verdict changed')

def positive_margin(decision,check=True):
    if decision['quartic_status']!='strictly_positive': raise ValueError('Need a strictly positive original input')
    if check:
        from verify_complete import audit_decision
        audit_decision(decision)
    eta=s.Rational(decision['barrier']['sample_epsilon'])
    return {'couplings':decision['couplings'],'quartic_status':'strictly_positive','normalized_lower_bound':str(5*eta/6),'componentwise_robustness_radius':str(5*eta/8),'inequality':'V4(K) > normalized_lower_bound * (Tr K)^2 for every nonzero physical configuration','radius_interpretation':'Every one of the 15 real coefficient components may change independently by strictly less than this radius; strict positivity is preserved.','complete_decision_audit_performed_here':check,'proof_dependency':'Requires correct complete positive classification of p-eta*q and Vq >=5/6 on the normalized physical domain; not a standalone positivity proof.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('decision');ap.add_argument('output');args=ap.parse_args();d=json.loads(Path(args.decision).read_text())
    if d['quartic_status']=='unstable':out=rational_witness(d)
    elif d['quartic_status']=='strictly_positive':out=positive_margin(d)
    else:raise ValueError('No nonzero robustness margin is asserted for an exactly flat potential')
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n');print(out.get('normalized_value',out.get('normalized_lower_bound')))
