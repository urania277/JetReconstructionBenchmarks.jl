# Julia Fastjet Utilities and Benchmarks

This repository has some utility and example functions for running jet
reconstruction with the
[`JetReconstruction`](https://github.com/JuliaHEP/JetReconstruction.jl) Julia
package.

## Overview

### `fastjet`

Example C++ application that runs fastjet for comparing with other
implementations:

- Timing comparisons
- Jet reconstruction outputs

Use the Julia script `Benchmark.jl` with the `--backend FastJet` option to run
systematic tests. Use the script `fastjet2json.jl` to convert fastjet outputs to
JSON files for comparison tests with Julia.

### `genevt`

Pythia8 application(s) for producing input HepMC3 files.

TODO: Proper CMake setup, but the source files will take you most of the way.

### `data`

Sorted HepMC3 data files used as reconstruction inputs (compressed).

Note that the `Benchmark.jl` script will trigger unpacking of the compressed
data files for running with fastjet.

### `src`

- `Benchmark.jl` run timing tests for different backends, allowing switching of
  algorithm, strategy, etc.
- `hepmc32summary.jl` write CSV files summarising the content of HepMC3 test
  input files (average cluster density)
- `fastjet2json.jl` convert fastjet outputs to JSON files for event validation
tests with Julia
- `generate-benchmarks-{pp,ee}.sh` example files of how to generate a set of
  benchmark files for various parameters
- `findmin-test.jl` tests the speed of finding minima in vectors, using Julia's
  builtin `findmin` and an turbo version that uses `LoopVectorization` (as is
  used in `JetReconstruction.jl`)

### `results`

Archive of results.

Recommended structure is to store in a particular subdirectory the results from
a particular machine. See, e.g., `generate-benchmarks.sh`.
