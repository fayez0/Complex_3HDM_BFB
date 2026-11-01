# StableWein comparison

The retained `stablewein_results.json` records an executed comparison against the public StableWein implementation using Wolfram Language 15.0.1.

Upstream source: `https://github.com/jurciukonis/StableWein`  
Pinned revision: `a540c6dbe099998f53d5af098e1c4a325c54ca84`

To rerun from this directory:

```bash
wolframscript -file run_stablewein.wls
```

The script downloads the pinned upstream notebook, checks its Wolfram-expression hash, evaluates the selected public Weinberg routines, and writes `stablewein_rerun_results.json`. Network access and a Wolfram kernel are required. The upstream notebook, binaries, and neural-network files are not redistributed here.

`verify_conventions.wls` checks the mapping between this repository's cyclic edge convention `(12,23,31)` and the StableWein convention. The executed comparison used the public modes 0--3 and the unmodified numerical minimizer on the six nonflat benchmarks. The five exactly flat examples are interpreted separately because StableWein uses strict quartic positivity, whereas this repository distinguishes nonnegative quartics with flat directions.

These are targeted validation records, not a representative accuracy measurement and not an independent proof of the universal criterion.
