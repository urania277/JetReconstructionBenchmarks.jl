#! /bin/sh
#
# Benchmark on different thread numbers

# Abort on error
set -e

max_threads=${1:-"10"}

echo "threads,time_per_event,event_rate"
threads=1
while [ $threads -le $max_threads ]; do
    julia --threads=$threads --project src/thread-run.jl -A AntiKt -R 0.4 --repeats 4 --nsamples 32 data/events-pp-13TeV-20GeV.hepmc3
    threads=$(($threads+1))
done
