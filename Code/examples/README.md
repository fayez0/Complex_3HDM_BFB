# Examples

The coefficient order is `(12,23,31)`. Each input has five arrays of length three:
`a`, `b`, `c`, `d_real`, and `d_imag`.

- `input_strictly_positive.json`: a nontrivial strictly positive regression input.
- `input_flat.json`: the identically zero quartic, used as a simple flat-boundary check.
- `input_unstable.json`: a complex-coefficient unstable benchmark.
- `Bminus_small_witness.json`: rational Higgs fields giving a strictly negative quartic value for the unstable benchmark; this can be checked with the Python standard library only.
- `positive_margin_example.json`: an example positive-margin certificate from the research validation.
- `unstable_witness_example.json`: an example instability-certificate record.

```bash
python weinberg_bfb.py classify examples/input_strictly_positive.json --output positive.decision.json
python weinberg_bfb.py classify examples/input_flat.json --output flat.decision.json
python weinberg_bfb.py classify examples/input_unstable.json --output unstable.decision.json
python -S weinberg_bfb.py verify-witness examples/Bminus_small_witness.json --output witness.audit.json
```
