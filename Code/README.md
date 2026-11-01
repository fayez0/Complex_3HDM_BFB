# Complex Weinberg 3HDM quartic stability

Exact research code accompanying the manuscript

> **A finite exact criterion for quartic stability in the complex Weinberg three-Higgs-doublet model**  
> Fayez Abu-Ajamieh (2026)

The code classifies the quartic potential of the complex `Z2 x Z2` Weinberg three-Higgs-doublet model as

- `strictly_positive`,
- `nonnegative_with_flat_direction`, or
- `unstable`.

Both neutral and charge-breaking field configurations are included. The implementation uses exact rational input components, algebraic root isolation, exact physicality tests, a fixed positive regularizer for exceptional coefficient strata, and a separate complete certificate audit.

## Potential and conventions

For `K_ij = Phi_i^dagger Phi_j` and `r_i = K_ii`,

```text
V4 = 1/2 sum_i a_i r_i^2
   + sum_(ij) [ b_ij r_i r_j + c_ij |K_ij|^2 + Re(d_ij K_ij^2) ].
```

The cyclic edge order is **(12,23,31)**. The physical normalized Gram domain is

```text
K >= 0,   Tr K = 1,   rank K <= 2.
```

The executable accepts five length-three arrays: `a`, `b`, `c`, `d_real`, `d_imag`. Use integers or exact rational/decimal **strings**; binary floating-point literals are deliberately rejected.

```json
{
  "a": [3, 4, 5],
  "b": [-1, 2, -2],
  "c": [2, 3, 1],
  "d_real": [1, 2, -1],
  "d_imag": [1, -1, 2]
}
```

## Installation

Python 3.13 and the exact-arithmetic dependencies in `requirements.txt` are sufficient for the public classifier:

```bash
python -m pip install -r requirements.txt
```

Run with assertions enabled; `python -O` and `python -OO` are intentionally rejected.

## Quick start

Classify an input and automatically run the complete audit:

```bash
python weinberg_bfb.py classify examples/input_strictly_positive.json --output decision.json
```

Re-audit a stored trusted decision:

```bash
python weinberg_bfb.py verify decision.json --output audit.json
```

Produce a physical certificate from an audited positive or unstable decision:

```bash
python weinberg_bfb.py physical decision.json --output physical_certificate.json
```

Check a rational negative-field witness using only the Python standard library:

```bash
python -S weinberg_bfb.py verify-witness examples/Bminus_small_witness.json --output witness_audit.json
```

A failed calculation, resource exhaustion, or unresolved algebraic step is **not** interpreted as a quartic-sign verdict.

## Scope

A successful `classify` output contains the three-way quartic status and a complete audit record. Instability can additionally be certified by explicit rational Higgs fields with `V4 < 0`; strict positivity can be accompanied by a certified lower margin.

The software classifies the **quartic part only**. If the quartic has flat directions, boundedness of the full quadratic-plus-quartic scalar potential is a separate question. The repository is for the symmetry-restricted complex Weinberg model, **not** the unrestricted general 3HDM.

## Repository structure

```text
weinberg_bfb.py          public command-line interface
src/                    exact classifier and independent auditors
tests/                  fast interface and exact-algebra regressions
examples/               representative inputs and physical certificates
validation/final/       analytically controlled reference checks
validation/stablewein/  retained external-comparison scripts and results
docs/CODE_MAP.md        map from the paper's construction to source modules
```

Early numerical prototypes, historical checkpoint archives, generated PDFs, and large intermediate decision files are intentionally omitted from this GitHub package.

## Tests and reproducibility

Fast regression suite:

```bash
python -m unittest discover -s tests -v
```

Independent exact identities for six analytic reference cases:

```bash
python validation/final/verify_analytic_references.py
```

Longer six-case exact classification driver:

```bash
python validation/final/run_closeout_cases.py
```

Some full exact classifications can take several minutes depending on the input and machine.

`validation/stablewein/` retains the scripts and records for the targeted comparison with the public StableWein implementation. StableWein uses strict quartic positivity, so exactly flat nonnegative potentials have a different convention from this repository's three-way output.

## Scientific status

This repository implements the finite criterion presented in the accompanying manuscript. The code and certificates have undergone internal exact-arithmetic, adversarial, cross-CAS, and targeted external comparisons. The universal completeness proof remains a mathematical claim in the paper rather than a formally verified theorem.

## Citation

Please cite the accompanying manuscript. `CITATION.cff` is included for GitHub's citation interface. Add the arXiv identifier and/or journal DOI once available.

## License

An open-source license has **not yet been selected**. See `LICENSE_STATUS.md` before making the repository public for reuse.
