#! /usr/bin/env julia
#
# Convert the text ouptout from the fastjet executables into a JSON
# format that is more suitable for integration into tests

using ArgParse
using Logging
using JSON

function read_fastjet_file(f)
    events = []
    jets = []
    event_number = 0
    for line in eachline(f)
        if startswith(line, "#") || !isnothing(match(r"^\s*$", line))
            continue
        end
        if startswith(line, "Jets in")
            # New event header
            if (event_number > 0) && (length(jets) > 0)
                push!(events, Dict("jetid" => event_number, "jets" => jets))
                jets = []
            end
            m = match(r"(\d+)", line)
            if !isnothing(m)
                event_number = parse(Int, m.match)
            else
                throw(ErrorException("Failed to find event number in line: $line"))
            end
        else
            # Jet output
            if !startswith(line, " ")
                continue
            end
            _, rap, phi, pt = split(line)
            push!(jets, Dict("rap" => parse(Float64, rap), "phi" => parse(Float64, phi), "pt" => parse(Float64, pt)))
        end
    end
    # Process last event
    push!(events, Dict("jetid" => event_number, "jets" => jets))
    events
end

parse_command_line(args) = begin
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
        "fastjetfile"
        help = "Fastjet output file to read"
        required = true

        "jsonoutput"
        help = "Output JSON filename"
        required = true
    end
    return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)

    event_dictonary = open(read_fastjet_file, args[:fastjetfile])
    JSON.print(open(args[:jsonoutput], write = true), event_dictonary, 2)
end

main()
