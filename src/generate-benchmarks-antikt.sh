#! /bin/sh
#
# Loop over parameters for generating benchmarks

# Abort on error
set -e

## Script parameters set from environment variables
# Generate output directory (adapt this to your setup)
output=${1:-"results/test-antikt"}
mkdir -p $output

echo "Benchmarking outputs: $output"

# Input datafile list
inputs=data/events-summary-ee-pp.csv

# Iterate over backends, including Python
for backend in JetReconstruction Fastjet AkTPython AkTNumPy; do
    for strategy in N2Tiled N2Plain; do
        for radius in 0.2 0.4 1.0 1.5 2.0 4.0; do
            # Python is slow, so do less trials
            if [ $backend == "JetReconstruction" -o $backend == "FastJet" ]; then
                trials=16
            else
                trials=4
            fi
            algorithm=AntiKt
            cmd="julia --project src/benchmark.jl --code $backend --radius $radius --algorithm $algorithm --strategy $strategy -m $trials --results $output $inputs"
            echo "Benchmark: $cmd"
            $cmd
        done
    done
done
