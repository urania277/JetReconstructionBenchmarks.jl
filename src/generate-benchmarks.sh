#! /bin/sh
#
# Loop over parameters for generating benchmarks

# Generate output directory (adapt this to your setup)
output=results/OSX-14.5-M2Pro-Julia-1.10.4
mkdir -p $output

# Input datafile list
inputs=data/events-summary-ee-pp.csv

# Iterate over backends
for backend in Julia FastJet; do
    for strategy in N2Plain N2Tiled; do
        for radius in 0.4 1.0; do
            echo Generating results for $backend at radius $radius using $strategy
            julia --project src/benchmark.jl --backend $backend -R $radius -A AntiKt -S $strategy -m 16 --results $output $inputs
        done
    done
done

