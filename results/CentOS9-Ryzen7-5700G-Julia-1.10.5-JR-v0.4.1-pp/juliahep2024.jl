### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 188be64e-3de6-11ef-222f-d7ff0fef6fe6
using DataFrames

# ╔═╡ 26572cc6-4e1e-47b3-9f9b-e32390d0798b
using CSV

# ╔═╡ 54b7f7ef-18c5-45ea-b8e8-acb68e63e27d
using Statistics

# ╔═╡ 451f22b7-c929-4c60-ab6c-0b3b174edbb6
using GLMakie

# ╔═╡ a5c2659e-ca5a-4b81-bead-6917a3146d86
using Logging

# ╔═╡ 24a8d5b4-c493-4965-aa66-983612d7d514
md"""# JuliaCON 2024 Plots

Plots for Jet Reconstruction talk at JuliaCON 2024
"""

# ╔═╡ 2a52599d-82b5-4632-b3ab-45e811611f8c
input_file = joinpath(@__DIR__, "..", "CentOS9-i7-Julia-1.10.4", "all-results.csv")

# ╔═╡ afba83db-a2e4-46ef-98de-4e1e6cde0d30
all_results_df = CSV.read(input_file, DataFrame);

# ╔═╡ d9cba1bf-c2f0-4305-9a84-3b3bb8feb5f1
all_results_df[1:5, :]

# ╔═╡ dc314a46-c44f-4941-b992-3a3a161e6b1e
names(all_results_df)

# ╔═╡ 4fc07c27-5d30-418a-bdec-b98f1efbd710
all_results_df[!, :backend] .== "Julia"

# ╔═╡ 3e09b13d-ea56-449d-b504-c9359494b5d0
names(all_results_df)

# ╔═╡ b2133112-6de4-44db-9c2f-1e523f26cee7
function to_array(arg)
    if !(typeof(arg) <: AbstractArray)
        return [arg]
    end
    arg
end

# ╔═╡ b0fcfe0d-6db4-4cfe-8b51-2a51a32ca9c6
"""
Primitive fixed selection criteria
"""
function select_results_rows(df::DataFrame; file = nothing, backend = nothing,
                             algorithm = nothing, strategy = nothing, R = nothing)
    selection = ones(Bool, size(all_results_df)[1])
    if !isnothing(file)
        file = to_array(file)
        selection = selection .& [x in file for x in all_results_df[!, :file]]
    end
    if !isnothing(backend)
        backend = to_array(backend)
        sel = [x in backend for x in all_results_df[!, :backend]]
        selection = selection .& sel
    end
    if !isnothing(algorithm)
        algorithm = to_array(algorithm)
        selection = selection .& [x in algorithm for x in all_results_df[!, :algorithm]]
    end
    if !isnothing(strategy)
        strategy = to_array(strategy)
        selection = selection .& [x in strategy for x in all_results_df[!, :strategy]]
    end
    if !isnothing(R)
        R = to_array(R)
        selection = selection .& [x in R for x in all_results_df[!, :R]]
    end
    selection
end

# ╔═╡ 3148ff99-50c0-4a47-b21e-93665b755877
Symbol("file")

# ╔═╡ b2a04431-a8fb-4cc3-a772-9bce310d1558
"""
    select_results_rows(df::DataFrame, selector::Dict)

Alternative selector that uses a `Dict` of selections for a lot
more flexibility
"""
function select_results_rows(df::DataFrame, selector::Dict)
    selection = ones(Bool, size(all_results_df)[1])
    for (col, value) in selector
        if !(col in names(df))
            @warn "$col is not a column name - ignored"
            continue
        end
        value = typeof(value) <: AbstractArray ? value : [value]
        selection = selection .& [x in value for x in df[!, col]]
    end
    selection
end

# ╔═╡ 3d801c3b-0add-4618-b342-a839c53e1ae7
select_results_rows(all_results_df; file = "events-ee-Z.hepmc3.gz")

# ╔═╡ 990602e0-fb4b-4459-b499-c6861e337b36
select_results_rows(all_results_df,
                    Dict("file" => "events-ee-Z.hepmc3.gz", "pants" => "wombat",
                         "R" => [0.2, 0.4]))

# ╔═╡ f7fc0d90-94f3-456b-9816-3d0a190f4ce7
md"## Plots"

# ╔═╡ aa24b6cd-a557-4ea3-9a17-36ebf9aef064
md"Julia Strategy vs Cluster Density"

# ╔═╡ f27502a6-98f9-42c1-9003-28c041940b57
julia_antikt_r04 = all_results_df[select_results_rows(all_results_df,
                                                      Dict("backend" => "Julia",
                                                           "algorithm" => "AntiKt",
                                                           "R" => 0.4)), :];

# ╔═╡ 41ce47d9-e324-4810-856d-893d7ce3bf6b
sort!(julia_antikt_r04, [:strategy, :R, :mean_particles]);

# ╔═╡ 01e58d51-02ea-412b-983d-aa99ecda2ff4
julia_antikt_r04_gdf = groupby(julia_antikt_r04, :strategy);

# ╔═╡ f6557478-abf4-4f93-ab48-a0d4805db752
begin
    f_julia_antikt_r04 = Figure()
    ax = Axis(f_julia_antikt_r04[1, 1],
              title = "Strategies - Julia Implementation ($(julia_antikt_r04_gdf[1][1, :algorithm]), R=$(julia_antikt_r04_gdf[1][1, :R]))",
              xlabel = "Average Cluster Density", ylabel = "μs per event")
    for (sn, slc) in enumerate(julia_antikt_r04_gdf)
        lines!(ax, slc[!, :mean_particles], slc[!, :time_per_event],
               color = Makie.wong_colors()[sn])
        scatter!(ax, slc[!, :mean_particles], slc[!, :time_per_event],
                 color = Makie.wong_colors()[sn], label = "$(slc[1, :strategy])")
    end
    axislegend(position = :rb)
end

# ╔═╡ 5db558d8-2b00-4926-898a-77314b88c23c
f_julia_antikt_r04

# ╔═╡ 315f02d8-32f1-4dd4-8604-dec659aa2147
save("julia-antikt-04.png", f_julia_antikt_r04)

# ╔═╡ 164ba10d-108f-4400-92de-7b92f8eac136
julia_antikt_r04_gdf[1][!, :time_per_event] ./ julia_antikt_r04_gdf[2][!, :time_per_event]

# ╔═╡ 44576d65-5ff5-4957-b003-387998fc75d1
md"Julia vs. FastJet - N2Plain"

# ╔═╡ 54e746d2-619f-475b-ac12-d18658933568
jfj_antikt = all_results_df[select_results_rows(all_results_df; algorithm = "AntiKt",
                                                strategy = "N2Plain", R = [0.4, 1.0]), :];

# ╔═╡ b01bebf0-321a-489d-8c69-a3208d2d5aac
sort!(jfj_antikt, [:backend, :R, :mean_particles]);

# ╔═╡ 5c318fa1-c094-4e72-84be-8a4cc6826785
jfj_antikt_gdf = groupby(jfj_antikt, [:backend, :R]);

# ╔═╡ b61c64f2-1ba6-4603-84ae-be4a821d766d
begin
    f_jfj_antikt = Figure()
    ax_jfj = Axis(f_jfj_antikt[1, 1],
                  title = "Julia - FastJet $(jfj_antikt_gdf[1][1, :strategy]) ($(jfj_antikt_gdf[1][1, :algorithm]))",
                  xlabel = "Average Cluster Density", ylabel = "μs per event")
    for (sn, slc) in enumerate(jfj_antikt_gdf)
        lines!(ax_jfj, slc[!, :mean_particles], slc[!, :time_per_event],
               color = Makie.wong_colors()[sn])
        scatter!(ax_jfj, slc[!, :mean_particles], slc[!, :time_per_event],
                 color = Makie.wong_colors()[sn],
                 label = "$(slc[1, :backend]) R=$(slc[1, :R])")
    end
    axislegend(position = :lt)
end

# ╔═╡ 96d54d97-5de3-402d-aa4d-0d5be8ebde68
f_jfj_antikt

# ╔═╡ 7a4b1d1a-aab9-4d47-a21f-3648e92c16bc
save("julia-fastjet-n2plain-04-10.png", f_jfj_antikt)

# ╔═╡ 64904420-b7e5-4f49-86a2-a28c8812416c
jfj_antikt[jfj_antikt[!, :file] .== "events-ee-Z.hepmc3.gz", :]

# ╔═╡ c12302e2-f894-4b2f-a5dd-bccc2b1def6c
md"Julia vs. FastJet N2Tiled"

# ╔═╡ 2b393e76-9fab-4a71-ac38-65194872e6ec
jfj_antikt_tiled = all_results_df[select_results_rows(all_results_df; algorithm = "AntiKt",
                                                      strategy = "N2Tiled", R = 0.4), :];

# ╔═╡ 5cec9213-1537-4169-b25f-1dc67c9a4e44
sort!(jfj_antikt_tiled, [:R, :backend, :mean_particles]);

# ╔═╡ 73cebe87-a853-40e0-8967-9519f56aba31
jfj_antikt_tiled_gdf = groupby(jfj_antikt_tiled, [:backend, :R]);

# ╔═╡ 08e38d26-27dd-4c60-b423-1f00088f8be1
jfj_antikt_tiled[(jfj_antikt_tiled[!, :file] .== "events-ee-Z.hepmc3.gz"), :]

# ╔═╡ 183cf56b-96ad-4bdd-98ac-d622c44a3fc4
begin
    f_jfj_antikt_tiled = Figure()
    ax_jfj_tiled = Axis(f_jfj_antikt_tiled[1, 1],
                        title = "Julia - FastJet $(jfj_antikt_tiled_gdf[1][1, :strategy]) ($(jfj_antikt_tiled_gdf[1][1, :algorithm]))",
                        xlabel = "Average Cluster Density", ylabel = "μs per event")
    for (sn, slc) in enumerate(jfj_antikt_tiled_gdf)
        lines!(ax_jfj_tiled, slc[!, :mean_particles], slc[!, :time_per_event],
               color = Makie.wong_colors()[sn])
        scatter!(ax_jfj_tiled, slc[!, :mean_particles], slc[!, :time_per_event],
                 color = Makie.wong_colors()[sn],
                 label = "$(slc[1, :backend]) R=$(slc[1, :R])")
    end
    axislegend(position = :lt)
end

# ╔═╡ 0506ec33-b39b-4dfb-8d49-779ab657816e
f_jfj_antikt_tiled

# ╔═╡ ee75b47a-6233-49a7-84e8-c2afe0198c37
save("julia-fastjet-n2tiled-04.png", f_jfj_antikt_tiled)

# ╔═╡ b8400920-d6a4-4e62-b01f-66b480494775
jfj_antikt_tiled_gdf[1][!, :time_per_event] ./ jfj_antikt_tiled_gdf[2][!, :time_per_event]

# ╔═╡ Cell order:
# ╠═24a8d5b4-c493-4965-aa66-983612d7d514
# ╠═188be64e-3de6-11ef-222f-d7ff0fef6fe6
# ╠═26572cc6-4e1e-47b3-9f9b-e32390d0798b
# ╠═54b7f7ef-18c5-45ea-b8e8-acb68e63e27d
# ╠═451f22b7-c929-4c60-ab6c-0b3b174edbb6
# ╠═a5c2659e-ca5a-4b81-bead-6917a3146d86
# ╠═2a52599d-82b5-4632-b3ab-45e811611f8c
# ╠═afba83db-a2e4-46ef-98de-4e1e6cde0d30
# ╠═d9cba1bf-c2f0-4305-9a84-3b3bb8feb5f1
# ╠═dc314a46-c44f-4941-b992-3a3a161e6b1e
# ╠═4fc07c27-5d30-418a-bdec-b98f1efbd710
# ╠═3e09b13d-ea56-449d-b504-c9359494b5d0
# ╠═b2133112-6de4-44db-9c2f-1e523f26cee7
# ╠═b0fcfe0d-6db4-4cfe-8b51-2a51a32ca9c6
# ╠═3148ff99-50c0-4a47-b21e-93665b755877
# ╠═b2a04431-a8fb-4cc3-a772-9bce310d1558
# ╠═3d801c3b-0add-4618-b342-a839c53e1ae7
# ╠═990602e0-fb4b-4459-b499-c6861e337b36
# ╟─f7fc0d90-94f3-456b-9816-3d0a190f4ce7
# ╟─aa24b6cd-a557-4ea3-9a17-36ebf9aef064
# ╠═f27502a6-98f9-42c1-9003-28c041940b57
# ╠═41ce47d9-e324-4810-856d-893d7ce3bf6b
# ╠═01e58d51-02ea-412b-983d-aa99ecda2ff4
# ╠═f6557478-abf4-4f93-ab48-a0d4805db752
# ╠═5db558d8-2b00-4926-898a-77314b88c23c
# ╠═315f02d8-32f1-4dd4-8604-dec659aa2147
# ╠═164ba10d-108f-4400-92de-7b92f8eac136
# ╟─44576d65-5ff5-4957-b003-387998fc75d1
# ╠═54e746d2-619f-475b-ac12-d18658933568
# ╠═b01bebf0-321a-489d-8c69-a3208d2d5aac
# ╠═5c318fa1-c094-4e72-84be-8a4cc6826785
# ╠═b61c64f2-1ba6-4603-84ae-be4a821d766d
# ╠═96d54d97-5de3-402d-aa4d-0d5be8ebde68
# ╠═7a4b1d1a-aab9-4d47-a21f-3648e92c16bc
# ╠═64904420-b7e5-4f49-86a2-a28c8812416c
# ╟─c12302e2-f894-4b2f-a5dd-bccc2b1def6c
# ╠═2b393e76-9fab-4a71-ac38-65194872e6ec
# ╠═5cec9213-1537-4169-b25f-1dc67c9a4e44
# ╠═73cebe87-a853-40e0-8967-9519f56aba31
# ╠═08e38d26-27dd-4c60-b423-1f00088f8be1
# ╠═183cf56b-96ad-4bdd-98ac-d622c44a3fc4
# ╠═0506ec33-b39b-4dfb-8d49-779ab657816e
# ╠═ee75b47a-6233-49a7-84e8-c2afe0198c37
# ╠═b8400920-d6a4-4e62-b01f-66b480494775
