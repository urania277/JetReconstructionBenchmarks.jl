# FastJet Test Applications

This code compiles a few small FastJet applications that are used for
benchmarking and validation of alternative implementations of clustering
algorithms.

The code requires the fastjet libraries (<https://fastjet.fr/>) as well as those
from HepMC3 (<https://gitlab.cern.ch/hepmc/HepMC3>).

Depending on your system setup you make need to add the path to the HepMC3 CMake
setup files to `CMAKE_PREFIX_PATH`.

## Compilation

Configure and compile using CMake in the standard way, e.g.,

```sh
cmake -S . -B build
cmake --build build
```

## Applications

### `fastjet-inclusive`

`fastjet-inclusive` will run fastjet with the standard suite of pp algorithms,
optionally outputting a list of inclusive jets. This is the standard executable
used to benchmark fastjet.

```sh
Usage: ./fastjet-inclusive [-h] [-m max_events] [-n trials] [-s strategy] [-p power] [-d dump_file] <HepMC3_input_file>
 HepMC3_input_file: File with input events in HepMC3 format
 -m max_events: default is -1, which is all the events in the file
 -n trials: default is 8, which is the number of repeats to do
 -s strategy: valid values are 'Best' (default), 'N2Plain', 'N2Tiled'
 -p power: -1=antikt, 0=cambridge_achen, 1=inclusive kt
 -R size: R parameter, cone size (default = 0.4)
 -P pt_min: minimum pt for inclusive jet output (default = 5.0)
 -d dump_file: output jets are printed to here (use '-' for stdout)
 -h: print this message
 ```

### `fastjet-exclusive`

`fastjet-exclusive` is used to generate exclusive jet outputs for validation.

### `fastjet2json.jl`

`fastjet2json.jl` converts the text output from the fastjet applications into a
JSON format, more suitable for integration tests.

```sh
julia --project=.. fastjet2json.jl fastjet_input_file json_output_file
```
