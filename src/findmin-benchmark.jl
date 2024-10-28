#! /usr/bin/env julia
#
# Test various strategies for finding the minimum value in an
# array of random numbers.
#
# This is a proxy for the jet reconstruction problem, where we
# require each iteration to search for the lowest dij value.

using Chairmarks
using Random
using LoopVectorization
using Statistics
using Printf
using SIMD

function fast_findmin(dij::DenseVector{T}, n) where T
    best = 1
    @inbounds dij_min = dij[1]
    @turbo for here in 2:n
        dij_here = dij[here]
        newmin = dij_here < dij_min
        best = newmin ? here : best
        dij_min = newmin ? dij_here : dij_min
    end
    dij_min, best
end

function basic_findmin(dij::DenseVector{T}, n) where T
    best = 1
    @inbounds dij_min = dij[1]
    @inbounds @simd for here in 2:n
        dij_here = dij[here]
        newmin = dij_here < dij_min
        best = newmin ? here : best
        dij_min = newmin ? dij_here : dij_min
        # best = ifelse(newmin, here, best)
        # dij_min = ifelse(newmin, dij_here, dij_min)
    end
    dij_min, best
end

function julia_findmin(dij::DenseVector{T}, n) where T
    findmin(@view dij[1:n])
end

function naive_findmin(dij::DenseVector{T}, n) where T
    x = @fastmath foldl(min, @view dij[1:n])
    i = findfirst(==(x), dij)::Int
    x, i
end

function naive_findmin_reduce(dij::DenseVector{T}, n) where T
    x = @fastmath reduce(min, @view dij[1:n])
    i = findfirst(==(x), dij)::Int
    x, i
end

function fast_findmin_simd(dij::DenseVector{T}, n) where T
    laneIndices = SIMD.Vec{8, Int}((1, 2, 3, 4, 5, 6, 7, 8))
    minvals = SIMD.Vec{8, T}(Inf)
    min_indices = SIMD.Vec{8, Int}(0)

    n_batches, remainder = divrem(n, 8)
    lane = VecRange{8}(0)
    i = 1
    @inbounds @fastmath for _ in 1:n_batches
        dijs = dij[lane + i]
        predicate = dijs < minvals
        minvals = vifelse(predicate, dijs, minvals)
        min_indices = vifelse(predicate, laneIndices, min_indices)

        i += 8
        laneIndices += 8
    end

    min_value = SIMD.minimum(minvals)
    min_index = @inbounds min_value == minvals[1] ? min_indices[1] : min_value == minvals[2] ? min_indices[2] :
                min_value == minvals[3] ? min_indices[3] : min_value == minvals[4] ? min_indices[4] :
                min_value == minvals[5] ? min_indices[5] : min_value == minvals[6] ? min_indices[6] :
                min_value == minvals[7] ? min_indices[7] : min_indices[8]

    @inbounds @fastmath for _ in 1:remainder
        xi = dij[i]
        pred = dij[i] < min_value
        min_value= ifelse(pred, xi, min_value)
        min_index = ifelse(pred, i, min_index)
        i += 1
    end
    return min_value, min_index
end

function run_descent(v::DenseVector{Float64}, f::T; perturb = 0) where T
    # Ensure we do something with the calculation to prevent the
    # compiler from optimizing everything away!
    sum = 0.0
    for n in length(v):-1:2
        val, _ = f(v, n)
        sum += val
        # If one wants to perturb the array do it like this, which
        # is a proxy for changing values as the algorithm progresses.
        for _ in 1:min(perturb,n)
            v[rand(1:n)] = abs(rand())
        end
    end
    sum
end

function report(f, v)
    bm = @be run_descent(v, f; perturb = 5)
    print("$(String(Symbol(f))) min: ", minimum(bm).time * 1e6, " μs; ")
    println("mean: ", mean(bm).time * 1e6, " μs")
end

function main(ARGS)
    # Setup the random array
    v = rand(400)

    # Run the benchmark
    println("Running the benchmark for an array of length $(length(v))")

    report(fast_findmin, v)
    report(basic_findmin, v)
    report(julia_findmin, v)
    report(naive_findmin, v)
    report(naive_findmin_reduce, v)
    report(fast_findmin_simd, v)
end

main(ARGS)
