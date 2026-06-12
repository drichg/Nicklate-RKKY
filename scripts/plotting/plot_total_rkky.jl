using JLD2
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function default_total_rkky_file()
	for filename in ("U=0-8-total.jld2", "U=0-8.0-total.jld2")
		filepath = joinpath(PROJECT_ROOT, "data", "results", filename)
		if isfile(filepath)
			return filepath
		end
	end

	return joinpath(PROJECT_ROOT, "data", "results", "U=0-8-total.jld2")
end

function plot_total_rkky_from_jld2(
	filepath::AbstractString = default_total_rkky_file();
	output_dir::AbstractString = joinpath(PROJECT_ROOT, "results", "rkky"),
)
	mkpath(output_dir)
	rslt = load(filepath, "rslt")
	U_values = collect(rslt.U_values)

	plots = Dict{String, Any}()
	for (name, values) in (
		"J0" => rslt.J0,
		"J1" => rslt.J1,
		"J2" => rslt.J2,
		"J3" => rslt.J3,
	)
		vector_values = collect(values)
		length(vector_values) == length(U_values) ||
			throw(ArgumentError("$name length must match U_values length"))

		plt = plot(
			U_values,
			vector_values;
			xlabel = "U",
			ylabel = name,
			title = "$name",
			legend = false,
			linewidth = 2,
			marker = :circle,
			markersize = 3,
		)

		savefig(plt, joinpath(output_dir, "rkky_total_$(name)_U0-8.png"))
		plots[name] = plt
	end

	return plots
end

if abspath(PROGRAM_FILE) == @__FILE__
	plot_total_rkky_from_jld2()
end
