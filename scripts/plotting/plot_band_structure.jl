function high_symmetry_path(n_per_segment::Integer = 100)
	n_per_segment > 0 || throw(ArgumentError("n_per_segment must be positive"))

	points = [(0.0, 0.0), (pi, 0.0), (pi, pi), (0.0, 0.0)]
	labels = ["Γ", "X", "M", "Γ"]

	kpoints = Vector{Tuple{Float64, Float64}}()
	distances = Float64[]
	ticks = Float64[]
	total_distance = 0.0

	push!(kpoints, points[1])
	push!(distances, total_distance)
	push!(ticks, total_distance)

	for segment in 1:(length(points)-1)
		start = points[segment]
		stop = points[segment+1]
		dkx = stop[1] - start[1]
		dky = stop[2] - start[2]
		segment_length = hypot(dkx, dky)

		for step in 1:n_per_segment
			fraction = step / n_per_segment
			kx = start[1] + fraction * dkx
			ky = start[2] + fraction * dky
			push!(kpoints, (kx, ky))
			push!(distances, total_distance + fraction * segment_length)
		end

		total_distance += segment_length
		push!(ticks, total_distance)
	end

	return reduce(vcat, ([kx ky] for (kx, ky) in kpoints)), distances, ticks, labels
end

function band_energies(kpoints, U, mu)
	bands = zeros(Float64, size(kpoints, 1), 4)

	for i in axes(kpoints, 1)
		H = tight_binding(kpoints[i, 1], kpoints[i, 2], U, mu)
		bands[i, :] = eigvals(Hermitian(H))
	end

	return bands
end

function plot_band_structure(
	U,
	mu;
	n_per_segment::Integer = 100,
	savepath::Union{Nothing, String} = joinpath(
		NickelateRKKY.PROJECT_ROOT,
		"results",
		"susceptibility",
		"band_structure.png",
	),
)
	kpoints, distances, ticks, labels = high_symmetry_path(n_per_segment)
	bands = band_energies(kpoints, U, mu)

	plt = plot(
		xlabel = "",
		ylabel = "Energy",
		legend = false,
		xticks = (ticks, labels),
		xlims = (first(distances), last(distances)),
	)

	for band in axes(bands, 2)
		plot!(plt, distances, bands[:, band], color = :black, linewidth = 1.8)
	end

	for tick in ticks
		vline!(plt, [tick], color = :gray, linestyle = :dash, linewidth = 0.8)
	end

	hline!(plt, [0.0], color = :gray, linestyle = :dot, linewidth = 0.8)

	if savepath !== nothing
		savefig(plt, savepath)
	end

	return plt, distances, bands
end

using LinearAlgebra
using Plots

include(joinpath(@__DIR__, "..", "..", "src", "NickelateRKKY.jl"))
using .NickelateRKKY

