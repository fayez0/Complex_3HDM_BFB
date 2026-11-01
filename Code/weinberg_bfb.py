#!/usr/bin/env python3
"""Author-review release v1.1: audited quartic-sign decisions and field certificates.

Coefficient inputs are strict JSON. Algebraic decision files are trusted research
artifacts, not a safe interchange format for arbitrary untrusted expressions.
"""
from __future__ import annotations
import argparse, json, os, sys, tempfile
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent
KEYS=('a','b','c','d_real','d_imag')
VERSION='1.1-author-review'


def _unique_object(pairs: list[tuple[str,Any]]) -> dict[str,Any]:
    out={}
    for key,value in pairs:
        if key in out: raise ValueError(f'Duplicate JSON key: {key}')
        out[key]=value
    return out


def load_json(path: str | Path) -> dict[str,Any]:
    with Path(path).open(encoding='utf-8') as f: data=json.load(f,object_pairs_hook=_unique_object)
    if not isinstance(data,dict): raise ValueError('Expected a JSON object')
    return data


def validate_coefficients(data: dict[str,Any]) -> dict[str,list[str]]:
    """Parse literals without evaluating expressions or silently accepting floats."""
    if set(data)!=set(KEYS): raise ValueError(f'Expected exactly these five arrays: {", ".join(KEYS)}')
    out={}
    for key in KEYS:
        row=data[key]
        if not isinstance(row,list) or len(row)!=3: raise ValueError(f'{key} must be an array of length three')
        parsed=[]
        for i,value in enumerate(row):
            if isinstance(value,bool) or not isinstance(value,(int,str)): raise ValueError(f'{key}[{i}] must be an integer or an exact rational/decimal string; binary floats and booleans are rejected')
            try: q=Fraction(value)
            except (ValueError,ZeroDivisionError,TypeError) as exc: raise ValueError(f'Invalid exact rational literal in {key}[{i}]: {value!r}') from exc
            parsed.append(str(q))
        out[key]=parsed
    return out


def write_json(path: str | Path,data: dict[str,Any]) -> None:
    """Atomically replace the destination only with a complete JSON record."""
    dest=Path(path); dest.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile('w',encoding='utf-8',dir=dest.parent,prefix='.'+dest.name+'.',suffix='.tmp',delete=False) as f:
        temp=Path(f.name); json.dump(data,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    try: os.replace(temp,dest)
    finally:
        if temp.exists(): temp.unlink()


def audit_any(record: dict[str,Any]) -> dict[str,Any]:
    if record.get('format')=='positive-scale-normalized-v1':
        from scale_normalized_sign import audit
    else:
        from robust_sign import audit
    return audit(record)


def main(argv: list[str] | None=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',action='version',version=VERSION)
    sub=parser.add_subparsers(dest='command',required=True)
    for name,helptext in (
        ('classify','Classify coefficients and run the full certificate audit'),
        ('verify','Recompute the complete audit of a trusted decision record'),
        ('physical','Produce an original-field negative witness or a certified positive margin'),
        ('verify-witness','Check original-field negativity with standard-library rational arithmetic')):
        p=sub.add_parser(name,help=helptext);p.add_argument('input');p.add_argument('--output',required=True);p.add_argument('--force',action='store_true',help='Replace an existing output; a failure replaces it with a record containing no physical verdict')
    args=parser.parse_args(argv);dest=Path(args.output)
    if dest.exists() and not args.force:
        print(json.dumps({'completed':False,'error':'Output exists; choose a new path or use --force'}),file=sys.stderr);return 2
    if Path(args.input).resolve()==dest.resolve():
        print(json.dumps({'completed':False,'error':'Input and output paths must differ'}),file=sys.stderr);return 2
    try:
        if not __debug__: raise RuntimeError('Optimized Python (-O/-OO) is unsupported: run with assertions enabled')
        raw=load_json(args.input);sys.path.insert(0,str(ROOT/'src'))
        if args.command=='classify':
            data=validate_coefficients(raw)
            from scale_normalized_sign import classify
            result=classify(data);result['complete_audit']=audit_any(result)
            if result['complete_audit'].get('verified') is not True: raise RuntimeError('Complete audit did not verify')
            result['release_interface']={'version':VERSION,'completed':True,'complete_audit_run':True}
        elif args.command=='verify':
            result=audit_any(raw)
            if result.get('verified') is not True: raise RuntimeError('Complete audit did not verify')
        elif args.command=='physical':
            verified=audit_any(raw)
            if verified.get('verified') is not True: raise RuntimeError('Source decision audit did not verify')
            if raw.get('format')=='positive-scale-normalized-v1':
                from scale_normalized_sign import physical_certificate
                result=physical_certificate(raw,check=False)
            else:
                from physical_certificates import rational_witness,positive_margin
                status=raw.get('quartic_status')
                if status=='unstable': result=rational_witness(raw)
                elif status=='strictly_positive': result=positive_margin(raw,check=False)
                else: raise ValueError('No nonzero robustness margin or general original-field flat witness is claimed for a flat input')
            result['complete_decision_audit_performed_here']=True
        else:
            from verify_rational_witness import verify
            result=verify(raw)
        write_json(dest,result);print(json.dumps({'completed':True,'output':str(dest),'result':result.get('quartic_status',result.get('status','verified'))}));return 0
    except (ValueError,KeyError,TypeError,ZeroDivisionError,ArithmeticError,NotImplementedError,AssertionError,RuntimeError,MemoryError,OSError) as exc:
        failure={'completed':False,'release_version':VERSION,'command':args.command,'error_type':type(exc).__name__,'error':str(exc),'scope':'No physical quartic-sign verdict is returned by this failed command.'}
        try:write_json(dest,failure)
        except OSError:pass
        print(json.dumps(failure),file=sys.stderr);return 2 if isinstance(exc,(ValueError,KeyError,TypeError,ZeroDivisionError)) else 3

if __name__=='__main__':raise SystemExit(main())
