sunny_runtime_version() = string(Base.pkgversion(Sunny))

function sunny_bilayer_q_path(cryst; n_points::Integer = 301)
    n_points >= 4 || throw(ArgumentError("n_points must be at least 4"))
    vertices = [
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.0, 0.0],
    ]
    labels = ["Gamma", "X", "M", "Gamma"]
    sunny_path = q_space_path(cryst, vertices, n_points; labels)
    q_points = reduce(vcat, (permutedims(collect(q)) for q in sunny_path.qs))
    distances = zeros(Float64, n_points)
    for index in 2:n_points
        delta = sunny_path.qs[index] - sunny_path.qs[index - 1]
        distances[index] = distances[index - 1] + norm(cryst.recipvecs * delta)
    end
    tick_indices = sunny_path.xticks[1]
    return (;
        sunny_path,
        vertices,
        labels,
        q_points,
        distances,
        tick_indices,
        tick_distances = distances[tick_indices],
    )
end

function _parse_instability_q(message)
    match_result = match(r"q\s*=\s*\[([^\]]+)\]", message)
    match_result === nothing && return nothing
    values = try
        parse.(Float64, strip.(split(match_result.captures[1], ',')))
    catch
        return nothing
    end
    return length(values) == 3 ? values : nothing
end

function classify_sunny_exception(exception; stage = "unknown")
    message = sprint(showerror, exception)
    unstable = occursin("Instability at wavevector", message) ||
        (occursin("Not an energy-minimum", message) && occursin("unstable", message))
    return (;
        status = unstable ? "unstable_spin_wave" : "software_error",
        stage = String(stage),
        exception_type = string(typeof(exception)),
        message,
        q = _parse_instability_q(message),
    )
end

function select_sunny_attempt(attempts)
    isempty(attempts) && throw(ArgumentError("attempts must not be empty"))
    stable = findall(attempt -> attempt.status == "ok", attempts)
    candidates = isempty(stable) ? eachindex(attempts) : stable
    energies = [attempts[index].final_energy_eV for index in candidates]
    return candidates[argmin(energies)]
end

function classify_spin_wave(energies_eV; negative_tolerance = 1e-8)
    all(isfinite, energies_eV) || return (;
        status = "unstable_spin_wave",
        minimum_frequency_eV = minimum(energies_eV),
        negative_count = 0,
        q_index = nothing,
        branch_index = nothing,
        reason = "non-finite frequency",
    )
    unstable_indices = findall(<(-negative_tolerance), energies_eV)
    minimum_index = argmin(energies_eV)
    branch_index, q_index = Tuple(minimum_index)
    return (;
        status = isempty(unstable_indices) ? "ok" : "unstable_spin_wave",
        minimum_frequency_eV = energies_eV[minimum_index],
        negative_count = length(unstable_indices),
        q_index,
        branch_index,
        reason = isempty(unstable_indices) ? "" : "negative frequency",
    )
end

function compute_sunny_bilayer_spectrum(
    sys;
    n_q::Integer = 301,
    negative_tolerance = 1e-8,
)
    path = sunny_bilayer_q_path(sys.crystal; n_points = n_q)
    try
        measure = ssf_perp(sys)
        bands = intensities_bands(
            SpinWaveTheory(sys; measure),
            path.sunny_path,
        )
        energies_eV = Array(bands.disp) ./ SUNNY_ENERGY_SCALE
        weights = Array(bands.data)
        all(isfinite, weights) || error("Sunny returned non-finite spectral weights")
        diagnostic = classify_spin_wave(energies_eV; negative_tolerance)
        return (;
            status = diagnostic.status,
            diagnostic,
            measure = "ssf_perp",
            energy_unit = "eV",
            energies_eV,
            weights,
            q_points = path.q_points,
            distances = path.distances,
            tick_indices = path.tick_indices,
            tick_distances = path.tick_distances,
            labels = path.labels,
            exception = nothing,
        )
    catch exception
        diagnostic = classify_sunny_exception(exception; stage = "spectrum")
        diagnostic.status == "unstable_spin_wave" || rethrow()
        return (;
            status = diagnostic.status,
            diagnostic = nothing,
            measure = "ssf_perp",
            energy_unit = "eV",
            energies_eV = zeros(Float64, 0, n_q),
            weights = zeros(Float64, 0, n_q),
            q_points = path.q_points,
            distances = path.distances,
            tick_indices = path.tick_indices,
            tick_distances = path.tick_distances,
            labels = path.labels,
            exception = diagnostic,
        )
    end
end

function plot_sunny_bilayer_spectrum(
    spectrum;
    U,
    filepath::AbstractString,
)
    isempty(spectrum.energies_eV) &&
        throw(ArgumentError("cannot plot a spectrum that aborted before returning bands"))
    mkpath(dirname(filepath))
    n_modes, n_q = size(spectrum.energies_eV)
    x = repeat(spectrum.distances; outer = n_modes)
    energies = vec(permutedims(spectrum.energies_eV))
    weights = max.(vec(permutedims(spectrum.weights)), 0.0)
    scale = maximum(weights)
    marker_sizes = iszero(scale) ? fill(2.0, length(weights)) :
        1.5 .+ 5.0 .* sqrt.(weights ./ scale)
    status_label = spectrum.status == "ok" ? "" : " UNSTABLE"
    plot_object = scatter(
        x,
        energies;
        marker_z = weights,
        markersize = marker_sizes,
        markerstrokewidth = 0,
        color = :viridis,
        colorbar_title = "ssf_perp",
        xlabel = "",
        ylabel = "Energy (eV)",
        title = "Bilayer spin waves, U=$(U) eV$(status_label)",
        legend = false,
        xticks = (spectrum.tick_distances, spectrum.labels),
        xlims = extrema(spectrum.distances),
    )
    for tick in spectrum.tick_distances
        vline!(
            plot_object,
            [tick];
            color = :gray,
            linestyle = :dash,
            linewidth = 0.7,
            label = false,
        )
    end
    hline!(plot_object, [0.0]; color = :black, linewidth = 0.7, label = false)
    savefig(plot_object, filepath)
    return plot_object
end

function compute_sunny_bilayer_attempt(
    U,
    exchanges;
    dims = (4, 4, 1),
    seed,
    n_q = 301,
    negative_tolerance = 1e-8,
)
    minimized = minimize_sunny_bilayer(;
        exchanges...,
        dims,
        seed,
    )
    spectrum = compute_sunny_bilayer_spectrum(
        minimized.sys;
        n_q,
        negative_tolerance,
    )
    return (;
        U = Float64(U),
        status = spectrum.status,
        dims,
        seed = Int(seed),
        initial_energy_eV = minimized.initial_energy_eV,
        final_energy_eV = minimized.final_energy_eV,
        energy_per_site_eV = minimized.energy_per_site_eV,
        converged = minimized.converged,
        spins = minimized.spins,
        phase = minimized.phase,
        layer_correlation = minimized.layer_correlation,
        spectrum,
    )
end

function run_sunny_point_with_retries(
    U,
    exchanges;
    initial_seed,
    n_q = 301,
    negative_tolerance = 1e-8,
    attempt_calculator = compute_sunny_bilayer_attempt,
)
    attempts = Any[]
    calculate(dims, seed) = attempt_calculator(
        U,
        exchanges;
        dims,
        seed,
        n_q,
        negative_tolerance,
    )

    push!(attempts, calculate((4, 4, 1), initial_seed))
    if attempts[end].status != "ok"
        for offset in 1:4
            push!(attempts, calculate((4, 4, 1), initial_seed + offset))
        end
    end
    if !any(attempt -> attempt.status == "ok", attempts)
        for offset in 101:102
            push!(attempts, calculate((6, 6, 1), initial_seed + offset))
        end
    end
    if !any(attempt -> attempt.status == "ok", attempts)
        for offset in 201:202
            push!(attempts, calculate((8, 8, 1), initial_seed + offset))
        end
    end
    selected_attempt = select_sunny_attempt(attempts)
    selected = attempts[selected_attempt]
    return (;
        U = Float64(U),
        status = selected.status,
        selected_attempt,
        attempts,
    )
end

function sunny_scan_metadata(
    rkky_result;
    base_seed,
    n_q,
    negative_tolerance,
)
    source = (;
        metadata = rkky_result.metadata,
        U_values = copy(rkky_result.completed_U_values),
        J0 = copy(rkky_result.J0),
        J1 = copy(rkky_result.J1),
        J2 = copy(rkky_result.J2),
        J3 = copy(rkky_result.J3),
        J1p = copy(rkky_result.J1p),
    )
    return (;
        schema_version = 1,
        sunny_version = sunny_runtime_version(),
        public_energy_unit = "eV",
        sunny_internal_energy_unit = "meV",
        lattice_constants_angstrom = (3.85, 3.85, 20.0),
        positions = [[0.0, 0.0, 0.25], [0.0, 0.0, 0.75]],
        spin = 0.5,
        g = 2.0,
        initial_dims = (4, 4, 1),
        base_seed = Int(base_seed),
        n_q = Int(n_q),
        negative_tolerance = Float64(negative_tolerance),
        path_vertices = [
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.0, 0.0, 0.0],
        ],
        path_labels = ["Gamma", "X", "M", "Gamma"],
        bonds = _sunny_bilayer_bonds(),
        source,
    )
end

function save_sunny_result_atomic(filepath::AbstractString, result)
    mkpath(dirname(filepath))
    temporary = filepath * ".tmp"
    isfile(temporary) && rm(temporary; force = true)
    try
        jldsave(temporary; result)
        mv(temporary, filepath; force = true)
    finally
        isfile(temporary) && rm(temporary; force = true)
    end
    return filepath
end

load_sunny_result(filepath::AbstractString) = load(filepath, "result")

function _validate_sunny_input(rkky_result)
    fields = (:completed_U_values, :J0, :J1, :J2, :J3, :J1p)
    lengths = [length(getproperty(rkky_result, field)) for field in fields]
    all(==(first(lengths)), lengths) ||
        throw(ArgumentError("RKKY arrays have inconsistent lengths"))
    all(isfinite, rkky_result.completed_U_values) ||
        throw(ArgumentError("RKKY U values must be finite"))
    for field in fields[2:end]
        all(isfinite, getproperty(rkky_result, field)) ||
            throw(ArgumentError("RKKY exchange values must be finite"))
    end
    return nothing
end

function _selected_attempt(point)
    return point.attempts[point.selected_attempt]
end

function plot_sunny_point(point, filepath::AbstractString)
    selected = _selected_attempt(point)
    if !isempty(selected.spectrum.energies_eV)
        return plot_sunny_bilayer_spectrum(
            selected.spectrum;
            U = point.U,
            filepath,
        )
    end
    mkpath(dirname(filepath))
    phase = selected.phase
    plot_object = heatmap(
        phase.structure_factor';
        xlabel = "q_x grid index",
        ylabel = "q_y grid index",
        title = "U=$(point.U) eV: $(point.status), phase=$(phase.phase)",
        colorbar_title = "S(q)",
        aspect_ratio = :equal,
    )
    savefig(plot_object, filepath)
    return plot_object
end

function _sunny_image_path(image_dir, U)
    label = replace(string(round(U; digits = 1)), "." => "p", "-" => "m")
    return joinpath(image_dir, "sunny_U$(label).png")
end

function scan_sunny_bilayer(
    rkky_result;
    output_path::AbstractString,
    image_dir::AbstractString,
    base_seed = 20260612,
    n_q = 301,
    negative_tolerance = 1e-8,
    overwrite = false,
    point_calculator = run_sunny_point_with_retries,
    plot_writer = plot_sunny_point,
)
    _validate_sunny_input(rkky_result)
    metadata = sunny_scan_metadata(
        rkky_result;
        base_seed,
        n_q,
        negative_tolerance,
    )
    result = if isfile(output_path) && !overwrite
        loaded = load_sunny_result(output_path)
        loaded.metadata == metadata ||
            throw(ArgumentError("Sunny resume metadata mismatch"))
        loaded
    else
        (;
            metadata,
            completed_U_values = Float64[],
            points = Any[],
        )
    end
    result.completed_U_values ==
        rkky_result.completed_U_values[1:length(result.completed_U_values)] ||
        throw(ArgumentError("Sunny completed U values are not a source prefix"))

    mkpath(image_dir)
    for point in result.points
        image_path = _sunny_image_path(image_dir, point.U)
        isfile(image_path) || plot_writer(point, image_path)
    end

    start_index = length(result.completed_U_values) + 1
    for index in start_index:length(rkky_result.completed_U_values)
        U = rkky_result.completed_U_values[index]
        exchanges = (;
            J0 = rkky_result.J0[index],
            J1 = rkky_result.J1[index],
            J2 = rkky_result.J2[index],
            J3 = rkky_result.J3[index],
            J1p = rkky_result.J1p[index],
        )
        point = point_calculator(
            U,
            exchanges;
            initial_seed = base_seed + index,
            n_q,
            negative_tolerance,
        )
        image_path = _sunny_image_path(image_dir, U)
        plot_writer(point, image_path)
        result = (;
            metadata,
            completed_U_values = [result.completed_U_values; U],
            points = [result.points; point],
        )
        save_sunny_result_atomic(output_path, result)
    end
    return result
end
