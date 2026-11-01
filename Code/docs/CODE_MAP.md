# Code map

The public entry point is `weinberg_bfb.py`.

- `scale_normalized_sign.py`: top-level exact classification with positive scale normalization.
- `robust_sign.py`: three-way sign classification and complete audit dispatch.
- `universal_sign.py`: finite-sample sign decision after certified deformation.
- `universal_regularizer.py`: fixed positive regularizer and regularity obstruction factors.
- `zero_barrier.py`: zero-energy obstruction and certified deformation radius.
- `separating_current.py`: exact neutral/charged full-support phase enumeration with separating projections.
- `exact_current.py`, `exact_phase_solver.py`: exact current equations, algebraic roots, and sign isolation.
- `verify_complete.py`: exhaustive certificate audit, including candidate/root coverage.
- `verify_sample.py`, `verify_barrier.py`: independent sample and barrier checks.
- `physical_certificates.py`: rational instability witnesses and positive margins.
- `verify_rational_witness.py`: standard-library-only check of a rational negative-field witness.
- `homogeneous_anchor_path.py`: exact homogeneous-anchor path identities.

Early numerical prototypes and historical checkpoints are intentionally omitted because they are not required by the final public classifier.
