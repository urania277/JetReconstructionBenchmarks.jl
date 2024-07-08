#! /usr/bin/env julia
#
# Take the many CSV files output by paramater scans of different jet reconstruction
# methods and merge them into a single CSV file.
#
using ArgParse
using CSV
using DataFrames
using Logging

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
        "inputs"
        help = "Input CSV files"
        nargs = '+'
        required = true

        "--output"
        help = "Output CSV file"
        required = true
    end
    return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)

    results_dictionary = Dict{AbstractString, DataFrame}()
    for input in args[:inputs]
        @info "Reading $input"
        df = CSV.read(input, DataFrame)
        results_dictionary[input] = df
    end

    # Now we setup a DataFrame which will hold both the numbers read from each file,
    # but also have columns for backend, algorithm, strategy and R.
    results_df=DataFrame(file=String[], 
        mean_particles=Float64[], 
        n_samples=Int[],
        time_per_event=Float64[],
        backend=String[],
        algorithm=String[],
        strategy=String[],
        R=Float64[]
        )
    metadata_re = r"^(\w+)-(\w+)-(\w+)-([\d\.]+)\.csv"
    for (input_file_path, one_result_df) in results_dictionary
        filename = basename(input_file_path)
        metadata_match = match(metadata_re, filename)
        if isnothing(metadata_match)
            @warn "Failed to match expected metadata for filename $filename"
            continue
        end
        backend, algorithm, strategy, R_str = metadata_match.captures
        R = tryparse(Float64, R_str)
        for row in eachrow(one_result_df)
            push!(results_df, (row[:File], 
                row[:mean_particles], 
                row[:n_samples],
                row[:time_per_event],
                backend,
                algorithm,
                strategy,
                R))
        end
    end
    @info "Writing merged outputs to $(args[:output])"
    CSV.write(args[:output], results_df)
end

main()
