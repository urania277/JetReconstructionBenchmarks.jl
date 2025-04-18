#! /usr/bin/env julia
#
# Take the many CSV files output by parameter scans of different jet
# reconstruction methods and merge them into a single CSV file.
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
        default = "all-results.csv"
    end
    return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)

    results_df = DataFrame()
    # We will read all the input files to DataFrame and append their
    # contents to the final merged DataFrame.
    for input in args[:inputs]
        if input == args[:output]
            @warn "Skipping re-reading of output file $(args[:output])"
            continue
        end
        @info "Reading $input"
        df = CSV.read(input, DataFrame)
        append!(results_df, df; cols = :intersect, promote=true)
    end

    @info "Writing merged outputs to $(args[:output])"
    CSV.write(args[:output], results_df)
end

main()
