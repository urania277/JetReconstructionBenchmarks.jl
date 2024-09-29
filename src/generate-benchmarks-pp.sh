#! /bin/sh
#
# Loop over parameters for generating benchmarks

# Abort on error
set -e

# Generate output directory (adapt this to your setup
# or set the JRB_OUTPUT variable)
output=results/${JRB_OUTPUT:-"test-run"}
mkdir -p $output

echo "Benchmarking outputs: $output"

# Input datafile list
inputs=data/events-summary-ee-pp.csv

# Iterate over backends
for backend in Julia FastJet; do
    for strategy in N2Plain N2Tiled; do
        for radius in 0.2 0.4 1.0 1.5 2.0 3.0; do
            for algorithm in AntiKt CA Kt; do
                if [ $algorithm == "AntiKt" ]; then
                    p=-1.0
                elif [ $algorithm == "CA" ]; then
                    p=0.0
                elif [ $algorithm == "Kt" ]; then
                    p=1.0
                fi
                cmd="julia --project src/benchmark.jl --backend $backend -p $p -R $radius -A $algorithm -S $strategy -m 16 --results $output $inputs"
                echo "Benchmark: $cmd"
                $cmd
            done
        done
    done
done
