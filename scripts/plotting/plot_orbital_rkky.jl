using JLD2
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function plot_rkky_from_jld2(
	filepath::AbstractString = joinpath(PROJECT_ROOT, "data", "results", "U=0-8.0.jld2");
	output_dir::AbstractString = joinpath(PROJECT_ROOT, "results", "rkky"),
)
	mkpath(output_dir)
	rslt = load(filepath, "rslt")
	U_values = rslt.U_values
	series_labels = ["x2", "z2", "x2 + z2"]

	function plot_series(name, values; savepath)
		size(values, 2) == length(series_labels) ||
			throw(ArgumentError("$name must have three columns"))
		size(values, 1) == length(U_values) ||
			throw(ArgumentError("$name length must match U_values length"))

		plt = plot(
			xlabel = "U",
			ylabel = name,
			title = name,
			legend = :best,
		)

		for col in axes(values, 2)
			plot!(
				plt,
				U_values,
				values[:, col];
				label = series_labels[col],
				line = (:solid, 2),
				marker = (:circle, 3),
				markerstrokewidth = 0,
			)
		end

		savefig(plt, savepath)
		return plt
	end

	plots = Dict{String, Any}()
	for (name, values) in (
		"J0" => rslt.J0,
		"J1" => rslt.J1,
		"J2" => rslt.J2,
		"J3" => rslt.J3,
	)
		plt = plot_series(name, values; savepath = joinpath(output_dir, "rkky_$(name)_U0-8.png"))
		plots[name] = plt
	end

	J3_over_J1 = rslt.J3 ./ rslt.J1
	plots["J3/J1"] = plot_series(
		"J3/J1",
		J3_over_J1;
		savepath = joinpath(output_dir, "rkky_J3_over_J1_U0-8.png"),
	)

	return plots,rslt
end

if abspath(PROGRAM_FILE) == @__FILE__
	plots, rslt = plot_rkky_from_jld2()
end
