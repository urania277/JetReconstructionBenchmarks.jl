using ArgParse
using JetReconstruction
using CSV
using DataFrames
# using Infiltrator

function makeJets(data::Vector{PseudoJet})
    finaljets = jet_reconstruct(data; p=0, R = 0.4)
    
    finaljets
end

"""
Read a file, reconstruct and select exclusive jets
"""
function get_exclusive_jet_selection(filename, ptmin, R)
    exclusive_jets = Vector{PseudoJet}[]
    cluster_seqs = ClusterSequence[]
    
    events = read_final_state_particles(filename)
    for event in events
        cs = jet_reconstruct(event; p=0, R = R)
        push!(exclusive_jets, inclusive_jets(cs; ptmin=ptmin, T=PseudoJet))
        push!(cluster_seqs, cs)
    end
    (cluster_seqs, exclusive_jets)
end

"""
Read all files in turn, constructing a dictionary of inclusive jets
"""
function process_input_files(filenames, ptmin, R)
    final_dict = Dict{String, Tuple{Vector{ClusterSequence}, Vector{Vector{PseudoJet}}}}()
    for file in filenames
        final_dict[file] = get_exclusive_jet_selection(file, ptmin, R)
    end
    final_dict
end

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
    "--filter"
    help = "Set filtering parameters"
    arg_type = JetFilter
    default = JetFilter(0.3, 3)
    
    "--trim"
    help = "Set trimming parameters"
    arg_type = JetTrim
    default = JetTrim(0.3, 0.3, JetAlgorithm.CA)
    
    "--massdrop"
    help = "Set mass drop parameters"
    arg_type = MassDropTagger
    default = MassDropTagger(0.67, 0.09)
    
    "--softdrop"
    help = "Set soft drop parameters"
    arg_type = SoftDropTagger
    default = SoftDropTagger(0.1, 2.0)
    
    "--ptmin"
    help = "pT cut for exclusive jet selection"
    arg_type = Float64
    default = 2.0
    
    "--distance", "-R"
    help = "Distance parameter for jet merging"
    arg_type = Float64
    default = 0.4
    
    "--nsamples", "-n"
    help = "Number of measurement points to acquire"
    arg_type = Int
    default = 16
    
    "--results"
    help = "Write results in CSV format to this file"
    
    "files"
    help = "Event files to process"
    required = true
    nargs = '+'
    
end
return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)
    
    if endswith(args[:files][1], ".csv")
        if length(args[:files]) > 1
            println("When CSV input is given, only one file can be used")
            exit(1)
        end
        hepmc3_files_df = CSV.read(args[:files][1], DataFrame)
        input_file_path = dirname(args[:files][1])
        input_files = hepmc3_files_df[:, :File]
        input_files_full = map(fname -> joinpath(input_file_path, fname), input_files)
        hepmc3_files_df[:, :File_path] = input_files_full
    else
        input_files_full = args[:files]
        hepmc3_files_df = DataFrame("File_path" => input_files_full)
        hepmc3_files_df[:, :File] = map(basename, hepmc3_files_df[:, :File_path])
    end
    
    excl_jets_dict = process_input_files(input_files_full, args[:ptmin], args[:distance])
    
    total_ejets_per_sample = Int[]
    for k in input_files_full
        jet_count = 0
        for event_jets in excl_jets_dict[k][2]
            # println(event_jets)
            jet_count += length(event_jets)
        end
        push!(total_ejets_per_sample, jet_count)
        @info "File $k, got $jet_count exclusive jets"
    end
    hepmc3_files_df[:, :ExJets] = total_ejets_per_sample
    
    mass_drop_sample = Float64[]
    soft_drop_sample = Float64[]
    filter_sample = Float64[]
    trim_sample = Float64[]

    for (f, (event_cseqs, event_jets)) in excl_jets_dict
        mass_drop_times = Float64[]
        soft_drop_times = Float64[]
        filtering = Float64[]
        trimming = Float64[]
        # @infiltrate
        n_ex_jets_v = hepmc3_files_df[hepmc3_files_df[!, :File_path] .== f, :ExJets]
        @assert length(n_ex_jets_v) == 1
        n_ex_jets = n_ex_jets_v[1]

        for _ in 1:args[:nsamples]
            md = sd = fi = tr = 0.0
            for (cs, ejets) in zip(event_cseqs, event_jets)
                fi += @elapsed [jet_filtering(jet, cs, args[:filter]) for jet in ejets]
                tr += @elapsed [jet_trimming(jet, cs, args[:trim]) for jet in ejets]
                md += @elapsed [mass_drop(jet, cs, args[:massdrop]) for jet in ejets]
                sd += @elapsed [soft_drop(jet, cs, args[:softdrop]) for jet in ejets]
            end
            push!(mass_drop_times, md * 10^6 / n_ex_jets)
            push!(soft_drop_times, sd * 10^6 / n_ex_jets)
            push!(filtering, fi * 10^6 / n_ex_jets)
            push!(trimming, tr * 10^6 / n_ex_jets)
        end
        push!(mass_drop_sample, minimum(mass_drop_times))
        push!(soft_drop_sample, minimum(soft_drop_times))
        push!(filter_sample, minimum(filtering))
        push!(trim_sample, minimum(trimming))
    end
    hepmc3_files_df[!, :filter_time] = filter_sample
    hepmc3_files_df[!, :trimming_time] = trim_sample
    hepmc3_files_df[!, :massdrop_time] = mass_drop_sample
    hepmc3_files_df[!, :softdrop_time] = soft_drop_sample
    
    if !isnothing(args[:results])
        CSV.write(args[:results], hepmc3_files_df)
    end
    
end


# filter_timing = Float64[]
# trimming_timing = Float64[]
# massdrop_timing = Float64[]
# softdrop_timing = Float64[]

# n_samples = Int[]

# for event_file in hepmc3_files_df[:, :File_path]
#     allEvents = read_final_state_particles(event_file)
#     event_num = length(allEvents)
#     avg_filter_time = 0.0
#     avg_trimming_time = 0.0
#     avg_massdrop_time = 0.0
#     avg_softdrop_time = 0.0
#     total_inclusive_jets = 0
#     warm_up = 0.0

#     push!(n_samples, args[:nsamples])

#     if args[:nsamples] > 1
#         @info "Doing initial warm-up run"
#         for event in allEvents
#             cluster = makeJets(event)
#             allJets = inclusive_jets(cluster; ptmin = 2.0, T = PseudoJet)
#             total_inclusive_jets += length(allJets)

#             warm_up += @elapsed [jet_filtering(jet, cluster, args[:filter]) for jet in allJets]
#             warm_up += @elapsed [jet_trimming(jet, cluster, args[:trim]) for jet in allJets]
#             warm_up += @elapsed [mass_drop(jet, cluster, args[:massdrop]) for jet in allJets]
#             warm_up += @elapsed [soft_drop(jet, cluster, args[:softdrop]) for jet in allJets]
#         end
#     end
#     println("File $event_file - got $total_inclusive_jets inclusive jets")
#     println("Warm up was $warm_up")

#     for j in 1:args[:nsamples]
#         global counts = 0
#         for event in allEvents
#             counts += length(event)
#             cluster = makeJets(event)
#             allJets = inclusive_jets(cluster; ptmin = 2.0, T = PseudoJet)

#             avg_filter_time += @elapsed [jet_filtering(jet, cluster, args[:filter]) for jet in allJets]
#             avg_trimming_time += @elapsed [jet_trimming(jet, cluster, args[:trim]) for jet in allJets]
#             avg_massdrop_time += @elapsed [mass_drop(jet, cluster, args[:massdrop]) for jet in allJets]
#             avg_softdrop_time += @elapsed [soft_drop(jet, cluster, args[:softdrop]) for jet in allJets]
#         end

#         counts /= event_num 
#     end

#     avg_filter_time /= args[:nsamples] * 10 ^ -4 
#     avg_trimming_time /= args[:nsamples] * 10 ^ -4
#     avg_massdrop_time /= args[:nsamples] * 10 ^ -4
#     avg_softdrop_time /= args[:nsamples] * 10 ^ -4

#     push!(filter_timing, avg_filter_time)
#     push!(trimming_timing, avg_trimming_time)
#     push!(massdrop_timing, avg_massdrop_time)
#     push!(softdrop_timing, avg_softdrop_time)

# end

# hepmc3_files_df[:, :n_samples] = n_samples
# hepmc3_files_df[:, :filter_time] = filter_timing
# hepmc3_files_df[:, :trimming_time] = trimming_timing
# hepmc3_files_df[:, :massdrop_time] = massdrop_timing
# hepmc3_files_df[:, :softdrop_time] = softdrop_timing

# println(hepmc3_files_df)

# if !isnothing(args[:results])
# 	CSV.write(args[:results], hepmc3_files_df)
# end
# end

main()
