using LinearAlgebra
using Plots

include(joinpath(@__DIR__, "..", "..", "src", "NickelateRKKY.jl"))
using .NickelateRKKY


function normalize_orbital(orbital)
	name = lowercase(strip(String(orbital)))

	if name == "x2" || name == "z2"
		return name
	end

	throw(ArgumentError("orbital must be \"x2\" or \"z2\""))
end

function orbital_chi_value(chi, orbital)
	name = normalize_orbital(orbital)

	if name == "x2"
		# return real(chi[1, 1] + chi[1, 3])
		return real(chi[1, 1])
	end

	# return real(chi[2, 2] + chi[2, 4])
	return real(chi[2, 2])
end

function orbital_label(orbital)
	name = normalize_orbital(orbital)
	return name == "x2" ? "chi11 + chi13" : "chi22 + chi24"
end

function high_symmetry_q_indices(Nk::Integer)
	Nk > 0 || throw(ArgumentError("Nk must be positive"))
	Nk % 4 == 0 || throw(ArgumentError("Nk must be a multiple of 4 so pi/2 and pi lie on the q grid"))

	half = div(Nk, 2)
	points = [(0, 0), (0, half), (half, half), (0, 0)]
	qindices = Vector{Tuple{Int, Int}}()

	push!(qindices, points[1])
	for segment in 1:(length(points)-1)
		start = points[segment]
		stop = points[segment+1]
		dix = stop[1] - start[1]
		diy = stop[2] - start[2]

		for step in 1:half
			x = start[1] + div(dix * step, half)
			y = start[2] + div(diy * step, half)
			push!(qindices, (x, y))
		end
	end

	return qindices
end

function high_symmetry_path(Nk::Integer)
	qindices = high_symmetry_q_indices(Nk)
	qpoints = Vector{Tuple{Float64, Float64}}()
	distances = Float64[]
	ticks = Float64[]
	labels = ["Gamma", "Y", "M", "P", "Gamma"]

	half = div(Nk, 2)
	quarter = div(Nk, 4)
	tick_indices = [(0, 0), (0, half), (half, half), (quarter, quarter)]

	total_distance = 0.0
	push!(qpoints, (0.0, 0.0))
	push!(distances, total_distance)
	push!(ticks, total_distance)

	for i in 2:length(qindices)
		prev = qindices[i-1]
		curr = qindices[i]
		qx = 2pi * curr[1] / Nk
		qy = 2pi * curr[2] / Nk
		prev_qx = 2pi * prev[1] / Nk
		prev_qy = 2pi * prev[2] / Nk

		total_distance += hypot(qx - prev_qx, qy - prev_qy)
		push!(qpoints, (qx, qy))
		push!(distances, total_distance)

		if curr in tick_indices[2:end] || (curr == (0, 0) && i == length(qindices))
			push!(ticks, total_distance)
		end
	end

	return reduce(vcat, ([qx qy] for (qx, qy) in qpoints)), distances, ticks, labels
end

function compute_chi_path(
	U,
	orbital;
	Nk::Integer = 100,
	temperature = 0.001,
	n_target = 1.5,
)
	name = normalize_orbital(orbital)
	qindices = high_symmetry_q_indices(Nk)
	qpoints, distances, ticks, labels = high_symmetry_path(Nk)
	mu = find_chemical_potential(n_target, U; Nk = Nk, temperature = temperature)
	cache = diagonalized_kmesh(Nk, U)
	values = zeros(Float64, length(qindices))

	for (i, (iqx, iqy)) in pairs(qindices)
		chi = chi_matrix_q_index(iqx, iqy, cache, mu; temperature = temperature)
		values[i] = orbital_chi_value(chi, name)
	end

	return (; U, orbital = name, mu, cache, qindices, qpoints, distances, ticks, labels, values)
end

function plot_chi_path(
	U,
	orbital;
	Nk::Integer = 100,
	temperature = 0.001,
	n_target = 1.5,
	savepath = nothing,
)
	data = compute_chi_path(
		U,
		orbital;
		Nk = Nk,
		temperature = temperature,
		n_target = n_target,
	)

	label = orbital_label(data.orbital)
	if savepath === nothing
		savepath = joinpath(
			NickelateRKKY.PROJECT_ROOT,
			"results",
			"susceptibility",
			"chi_$(data.orbital)_U$(replace(string(U), "-" => "m")).png",
		)
	end

	plt = plot(
		data.distances,
		data.values;
		xlabel = "",
		ylabel = "Re [$(label)]",
		title = "Re [$(label)], U=$(U)",
		legend = false,
		xticks = (data.ticks, data.labels),
		xlims = (first(data.distances), last(data.distances)),
		# ylims = (0, 1.2),
		linewidth = 2,
	)

	for tick in data.ticks
		vline!(plt, [tick]; color = :gray, linestyle = :dash, linewidth = 0.8)
	end

	savefig(plt, savepath)
	println("mu = ", data.mu)
	println("saved ", savepath)

	return plt, data
end

function plot_chi_x2_multi_u(
	U_values;
	Nk::Integer = 100,
	temperature = 0.01,
	n_target = 1.5,
	offset = 0.2,
	savepath = joinpath(
		NickelateRKKY.PROJECT_ROOT,
		"results",
		"susceptibility",
		"chi_x2_multi_U.png",
	),
)
	U_list = collect(U_values)
	isempty(U_list) && throw(ArgumentError("U_values must not be empty"))
	isfinite(offset) || throw(ArgumentError("offset must be finite"))

	datasets = [
		compute_chi_path(
			U,
			"x2";
			Nk = Nk,
			temperature = temperature,
			n_target = n_target,
		) for U in U_list
	]

	first_data = first(datasets)
	plt = plot(
		xlabel = "",
		ylabel = "Re [chi11 + chi13] + offset",
		title = "x2 susceptibility for different U",
		xticks = (first_data.ticks, first_data.labels),
		xlims = (first(first_data.distances), last(first_data.distances)),
		linewidth = 2,
	)

	for (i, data) in pairs(datasets)
		shifted_values = data.values .+ (i - 1) * offset
		plot!(
			plt,
			data.distances,
			shifted_values;
			label = "U=$(data.U)",
			linewidth = 2,
		)
	end

	for tick in first_data.ticks
		vline!(plt, [tick]; color = :gray, linestyle = :dash, linewidth = 0.8, label = false)
	end

	if savepath !== nothing
		savefig(plt, savepath)
		println("saved ", savepath)
	end

	for data in datasets
		println("U = ", data.U, ", mu = ", data.mu)
	end

	return plt, datasets
end
function plot_chi_x2_multi_u_colormap(
	U_values;
	Nk::Integer = 100,
	temperature = 0.01,
	n_target = 1.5,
	offset = 0.2,
	color_scheme = :viridis,
	savepath = joinpath(
		NickelateRKKY.PROJECT_ROOT,
		"results",
		"susceptibility",
		"chi_x2_multi_U_colormap.png",
	),
)
	U_list = collect(U_values)
	isempty(U_list) && throw(ArgumentError("U_values must not be empty"))
	isfinite(offset) || throw(ArgumentError("offset must be finite"))

	datasets = [
		compute_chi_path(
			U,
			"x2";
			Nk = Nk,
			temperature = temperature,
			n_target = n_target,
		) for U in U_list
	]

	u_min, u_max = extrema(Float64.(U_list))
	u_span = u_max - u_min
	color_gradient = cgrad(color_scheme)
	first_data = first(datasets)

	plt = plot(
		xlabel = "",
		ylabel = "χ_x2 ",
		title = "x2 susceptibility for different U",
		xticks = (first_data.ticks, first_data.labels),
		xlims = (first(first_data.distances), last(first_data.distances)),
		legend = false,
	)

	for (i, data) in pairs(datasets)
		color_position = iszero(u_span) ? 0.5 : (Float64(data.U) - u_min) / u_span
		plot!(
			plt,
			data.distances,
			data.values .+ (i - 1) * offset;
			color = color_gradient[color_position],
			linewidth = 2,
			label = false,
		)
	end

	for tick in first_data.ticks
		vline!(plt, [tick]; color = :gray, linestyle = :dash, linewidth = 0.8, label = false)
	end

	scatter!(
		plt,
		[NaN, NaN],
		[NaN, NaN];
		marker_z = [u_min, u_max],
		color = color_scheme,
		clims = (u_min, u_max),
		colorbar = true,
		colorbar_title = "U",
		markersize = 0,
		label = false,
	)

	if savepath !== nothing
		savefig(plt, savepath)
		println("saved ", savepath)
	end

	for data in datasets
		println("U = ", data.U, ", mu = ", data.mu)
	end

	return plt, datasets
end
function plot_chi_z2_multi_u_colormap(
	U_values;
	Nk::Integer = 100,
	temperature = 0.01,
	n_target = 1.5,
	offset = 0.2,
	color_scheme = :viridis,
	savepath = joinpath(
		NickelateRKKY.PROJECT_ROOT,
		"results",
		"susceptibility",
		"chi_z2_multi_U_colormap.png",
	),
)
	U_list = collect(U_values)
	isempty(U_list) && throw(ArgumentError("U_values must not be empty"))
	isfinite(offset) || throw(ArgumentError("offset must be finite"))

	datasets = [
		compute_chi_path(
			U,
			"z2";
			Nk = Nk,
			temperature = temperature,
			n_target = n_target,
		) for U in U_list
	]

	u_min, u_max = extrema(Float64.(U_list))
	u_span = u_max - u_min
	color_gradient = cgrad(color_scheme)
	first_data = first(datasets)

	plt = plot(
		xlabel = "",
		ylabel = "χ_z2 ",
		title = "z2 susceptibility for different U",
		xticks = (first_data.ticks, first_data.labels),
		xlims = (first(first_data.distances), last(first_data.distances)),
		legend = false,
	)

	for (i, data) in pairs(datasets)
		color_position = iszero(u_span) ? 0.5 : (Float64(data.U) - u_min) / u_span
		plot!(
			plt,
			data.distances,
			data.values .+ (i - 1) * offset;
			color = color_gradient[color_position],
			linewidth = 2,
			label = false,
		)
	end

	for tick in first_data.ticks
		vline!(plt, [tick]; color = :gray, linestyle = :dash, linewidth = 0.8, label = false)
	end

	scatter!(
		plt,
		[NaN, NaN],
		[NaN, NaN];
		marker_z = [u_min, u_max],
		color = color_scheme,
		clims = (u_min, u_max),
		colorbar = true,
		colorbar_title = "U",
		markersize = 0,
		label = false,
	)

	if savepath !== nothing
		savefig(plt, savepath)
		println("saved ", savepath)
	end

	for data in datasets
		println("U = ", data.U, ", mu = ", data.mu)
	end

	return plt, datasets
end

function normalize_chi_channel(channel)
	name = lowercase(strip(String(channel)))
	valid_channels = ("11", "13", "x2_even", "x2_odd", "z2_even", "z2_odd")

	name in valid_channels ||
		throw(ArgumentError("channel must be one of: $(join(valid_channels, ", "))"))

	return name
end

function chi_channel_value(chi, channel)
	name = normalize_chi_channel(channel)

	if name == "11"
		return real(chi[1, 1])
	elseif name == "13"
		return real(chi[1, 3])
	elseif name == "x2_even"
		return real(chi[1, 1] + chi[1, 3])
	elseif name == "x2_odd"
		return real(chi[1, 1] - chi[1, 3])
	elseif name == "z2_even"
		return real(chi[2, 2] + chi[2, 4])
	end

	return real(chi[2, 2] - chi[2, 4])
end

function compute_chi_channel_path(
	U,
	channel;
	Nk::Integer = 100,
	temperature = 0.01,
	n_target = 1.5,
)
	name = normalize_chi_channel(channel)
	qindices = high_symmetry_q_indices(Nk)
	qpoints, distances, ticks, labels = high_symmetry_path(Nk)
	mu = find_chemical_potential(n_target, U; Nk = Nk, temperature = temperature)
	cache = diagonalized_kmesh(Nk, U)
	values = zeros(Float64, length(qindices))

	for (i, (iqx, iqy)) in pairs(qindices)
		chi = chi_matrix_q_index(iqx, iqy, cache, mu; temperature = temperature)
		values[i] = chi_channel_value(chi, name)
	end

	return (;
		U,
		channel = name,
		mu,
		cache,
		qindices,
		qpoints,
		distances,
		ticks,
		labels,
		values,
	)
end

function plot_chi_multi_u_colormap(
	U_values,
	channel;
	Nk::Integer = 100,
	temperature = 0.01,
	n_target = 1.5,
	offset = 0.2,
	color_scheme = :viridis,
	savepath = nothing,
)
	U_list = collect(U_values)
	isempty(U_list) && throw(ArgumentError("U_values must not be empty"))
	isfinite(offset) || throw(ArgumentError("offset must be finite"))
	name = normalize_chi_channel(channel)

	datasets = [
		compute_chi_channel_path(
			U,
			name;
			Nk = Nk,
			temperature = temperature,
			n_target = n_target,
		) for U in U_list
	]

	u_min, u_max = extrema(Float64.(U_list))
	u_span = u_max - u_min
	color_gradient = cgrad(color_scheme)
	first_data = first(datasets)

	plt = plot(
		xlabel = "",
		ylabel = "Re chi_$(name) + offset",
		title = "$(name) susceptibility for different U",
		xticks = (first_data.ticks, first_data.labels),
		xlims = (first(first_data.distances), last(first_data.distances)),
		legend = false,
	)

	for (i, data) in pairs(datasets)
		color_position = iszero(u_span) ? 0.5 : (Float64(data.U) - u_min) / u_span
		plot!(
			plt,
			data.distances,
			data.values .+ (i - 1) * offset;
			color = color_gradient[color_position],
			linewidth = 2,
			label = false,
		)
	end

	for tick in first_data.ticks
		vline!(plt, [tick]; color = :gray, linestyle = :dash, linewidth = 0.8, label = false)
	end

	scatter!(
		plt,
		[NaN, NaN],
		[NaN, NaN];
		marker_z = [u_min, u_max],
		color = color_scheme,
		clims = (u_min, u_max),
		colorbar = true,
		colorbar_title = "U",
		markersize = 0,
		label = false,
	)

	if savepath === nothing
		savepath = joinpath(
			NickelateRKKY.PROJECT_ROOT,
			"results",
			"susceptibility",
			"chi_$(name)_multi_U_colormap.png",
		)
	end
	savefig(plt, savepath)
	println("saved ", savepath)

	for data in datasets
		println("U = ", data.U, ", mu = ", data.mu)
	end

	return plt, datasets
end

function compute_chi_phase_heatmap(
	U,
	channel,
	R::Tuple{<:Integer, <:Integer};
	Nk::Integer = 100,
	temperature = 0.01,
	n_target = 1.5,
)
	Nk > 0 || throw(ArgumentError("Nk must be positive"))
	iseven(Nk) || throw(ArgumentError("Nk must be even for a centered momentum grid"))

	name = normalize_chi_channel(channel)
	mu = find_chemical_potential(n_target, U; Nk = Nk, temperature = temperature)
	cache = diagonalized_kmesh(Nk, U)
	q_indices = collect((-div(Nk, 2)):(div(Nk, 2)-1))
	q_values = 2pi .* q_indices ./ Nk
	chi_values = zeros(Float64, Nk, Nk)
	phase_factor = zeros(Float64, Nk, Nk)
	phase_values = zeros(Float64, Nk, Nk)

	for (ix, iqx) in pairs(q_indices), (iy, iqy) in pairs(q_indices)
		chi = chi_matrix_q_index(iqx, iqy, cache, mu; temperature = temperature)
		chi_value = chi_channel_value(chi, name)
		phase = cos(q_values[ix] * R[1] + q_values[iy] * R[2])
		chi_values[iy, ix] = chi_value
		phase_factor[iy, ix] = phase
		phase_values[iy, ix] = phase * chi_value
	end

	return (;
		U,
		channel = name,
		R,
		mu,
		cache,
		q_indices,
		q_values,
		chi_values,
		phase_factor,
		phase_values,
	)
end

function plot_chi_phase_heatmap(
	U,
	channel,
	R::Tuple{<:Integer, <:Integer};
	Nk::Integer = 100,
	temperature = 0.01,
	n_target = 1.5,
	color_scheme = :balance,
	savepath = nothing,
)
	data = compute_chi_phase_heatmap(
		U,
		channel,
		R;
		Nk = Nk,
		temperature = temperature,
		n_target = n_target,
	)

	if savepath === nothing
		savepath = joinpath(
			NickelateRKKY.PROJECT_ROOT,
			"results",
			"susceptibility",
			"chi_phase_$(data.channel)_R$(R[1])_$(R[2])_U$(replace(string(U), "-" => "m")).png",
		)
	end

	color_limit = maximum(abs, data.phase_values)
	color_limit = iszero(color_limit) ? 1.0 : color_limit
	ticks = ([-pi, 0.0, pi], ["-pi", "0", "pi"])
	plt = heatmap(
		data.q_values,
		data.q_values,
		data.phase_values;
		xlabel = "qx",
		ylabel = "qy",
		title = "cos(q.R) Re chi_$(data.channel)(q), R=$(R), U=$(U)",
		xticks = ticks,
		yticks = ticks,
		xlims = (-pi, pi),
		ylims = (-pi, pi),
		aspect_ratio = :equal,
		color = color_scheme,
		clims = (-color_limit, color_limit),
		colorbar_title = "cos(q.R) Re chi",
	)

	if savepath !== nothing
		savefig(plt, savepath)
		println("mu = ", data.mu)
		println("saved ", savepath)
	end

	return plt, data
end

# if abspath(PROGRAM_FILE) == @__FILE__
# Edit these values, then run:
# julia --project=. plot_chi_highsym.jl
# U_values = [0.0, 1.0, 2.0]
# Nk = 100
# temperature = 0.01
# n_target = 1.5
# offset = 0.1
# savepath = "chi_z2_multi_U.png"

# plt, datasets = plot_chi_x2_multi_u(
# 	U_values;
# 	Nk = Nk,
# 	temperature = temperature,
# 	n_target = n_target,
# 	offset = offset,
# 	savepath = savepath,
# )

# display(plt)
# end
# plt, datasets = plot_chi_z2_multi_u_colormap(
# 	0.0:0.5:8.0;
# 	Nk = 100,
# 	temperature = 0.01,
# 	offset = offset,
# 	color_scheme = :viridis,
# 	savepath = "chi_x2_colormap.png",
# )

# display(plt)
# plt, data = plot_chi_phase_heatmap(
# 	6.0,
# 	"x2_even",
# 	(2, 0);
# 	Nk = 100,
# 	temperature = 0.01,
# )

# display(plt)
if abspath(PROGRAM_FILE) == @__FILE__
	plot_chi_multi_u_colormap(
		0.0:0.5:8.0,
		"x2_odd";
		offset = 0.1,
		savepath = joinpath(
			NickelateRKKY.PROJECT_ROOT,
			"results",
			"susceptibility",
			"chi_x2_odd_multi_U_colormap.png",
		),
	)
end
