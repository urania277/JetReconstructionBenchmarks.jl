#! /usr/bin/env julia

using ArgParse
using LorentzVectorHEP
using Statistics
using UnicodePlots
using Logging
using CodecZlib

struct Particle{T}
    momentum::LorentzVector{T}
    status::Integer
    pdgid::Integer
    barcode::Integer
    vertex::Integer
end

struct FileSummary
    filename::String
    average_density::Float64
end

Particle{T}() where T = Particle(LorentzVector{T}(0., 0., 0., 0.), 0, 0, 0, 0)

""" Read a [HepMC3](https://doi.org/10.1016/j.cpc.2020.107310) ascii file.

    Each event is passed to the provided function f as a vector of Particles. A
    maximum number of events to read (value -1 to read all availble events) and
    a number of events to skip at the beginning of the file can be provided.
"""
function read_events(f, fin; maxevents=-1, skipevents=0)
    T = Float64
    particles = Particle{T}[]
    ievent = 0
    ipart = 0
    toskip = skipevents
    
    for (il, l) in enumerate(eachline(fin))
        if occursin(r"HepMC::.*-END_EVENT_LISTING", l)
            break
        end

        tok = split(l)

        if tok[1] == "E"
            ievent += 1
            (maxevents >= 0 && ievent > maxevents) && break
            if ievent > 1 && toskip == 0
                f(particles)
            end
            if toskip > 0
                toskip -= 1
            end
            resize!(particles, 0)
            ipart = 0
        elseif tok[1] == "P" && toskip == 0
            ipart += 1
            if ipart > length(particles)
                push!(particles, Particle{T}())
            end
            barcode = parse(Int, tok[2])
            vertex = parse(Int, tok[3])
            pdgid = parse(Int, tok[4])
            px = parse(T, tok[5])
            py = parse(T, tok[6])
            pz = parse(T, tok[7])
            e =  parse(T, tok[8])
            status = parse(Int, tok[10])
            push!(particles, Particle{T}(LorentzVector(e,px,py,pz), status, pdgid, barcode, vertex))
        end
    end
    #processing the last event:
    ievent > 0 && f(particles)
end


read_events(fname, maxevents=-1, skipevents=0) = begin
    f = open(fname)
    if endswith(fname, ".gz")
        @debug "Reading gzipped file $fname"
        f = GzipDecompressorStream(f)
    end

    events = Vector{LorentzVector}[]

    read_events(f, maxevents=maxevents, skipevents=skipevents) do parts
        input_particles = LorentzVector[]
        for p in parts
            if p.status == 1
                push!(input_particles, p.momentum)
            end
        end
        push!(events, input_particles)
    end

    events
end

parse_command_line(args) = begin
	s = ArgParseSettings(autofix_names = true)
	@add_arg_table! s begin
        "--summary"
        help = "Print only summary information, filename and average density"
        action = :store_true

		"--maxevents", "-n"
		help = "Maximum number of events to read. -1 to read all events from the  file."
		arg_type = Int
		default = -1

		"--skip", "-s"
		help = "Number of events to skip at beginning of the file."
		arg_type = Int
		default = 0

		"--info"
		help = "Print info level log messages"
		action = :store_true

		"--debug"
		help = "Print debug level log messages"
		action = :store_true

		"files"
		help = "The HepMC3 event files to read."
		required = true
        nargs = '+'
	end
	return parse_args(args, s; as_symbols = true)
end


function main()
    args = parse_command_line(ARGS)
    if args[:debug]
		logger = ConsoleLogger(stdout, Logging.Debug)
	elseif args[:info]
		logger = ConsoleLogger(stdout, Logging.Info)
	else
		logger = ConsoleLogger(stdout, Logging.Warn)
	end
	global_logger(logger)
    
    results = FileSummary[]

    for file in args[:files]
        events = read_events(file, args[:maxevents], args[:skip])
        n_events = length(events)
        n_particles = Int[]
        for (i, e) ∈ enumerate(events)
            push!(n_particles, length(e))
            @info "Event $(i) - $(length(e))"
        end
        average_n = mean(n_particles)
        push!(results, FileSummary(basename(file), average_n))
        if !args[:summary]
            println("File $file")
            println("  Number of events: $n_events")
            println("  Average number of particles: ", mean(n_particles))
            println(histogram(n_particles))
        end
    end

    println("File,mean_particles")
    sort!(results, by = x -> x.average_density)
    for r in results
        println(r.filename, ",", r.average_density)
    end
end

main()
