#! /bin/sh
#
# Loop over parameters for generating benchmarks

# Abort on error
set -e

# Generate output directory (adapt this to your setup)
output=results/${JRB_OUTPUT:-"test-run"}
mkdir -p $output

echo "Benchmarking outputs: $output"

# Input datafile list
inputs=data/events-summary-ee.csv

# Iterate over backends
for backend in Julia FastJet; do
    # Durham has no R or p parameter (or implicitly R=4, p=1)
    algorithm=Durham
    cmd="julia --project src/benchmark.jl --backend $backend -A $algorithm -m 16 --results $output $inputs"
    echo "Benchmark: $cmd"
    $cmd
    for radius in 0.2 0.4 1.0 1.5 2.0 4.0; do
        algorithm=EEKt
        for power in -1.0 0.0 1.0; do
            cmd="julia --project src/benchmark.jl --backend $backend -R $radius -A $algorithm -p $power -m 16 --results $output $inputs"
            echo "Benchmark: $cmd"
            $cmd
        done
    done
done
