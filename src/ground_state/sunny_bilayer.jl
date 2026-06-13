const SUNNY_ENERGY_SCALE = 1000.0

function _sunny_bilayer_bonds()
    return (;
        J0 = [(1, 2, (0, 0, 0))],
        J1 = [(1, 1, (1, 0, 0)), (2, 2, (1, 0, 0))],
        J2 = [(1, 1, (1, 1, 0)), (2, 2, (1, 1, 0))],
        J3 = [(1, 1, (2, 0, 0)), (2, 2, (2, 0, 0))],
        J1p = [(1, 2, (1, 0, 0))],
    )
end

function build_sunny_bilayer(;
    J0,
    J1,
    J2,
    J3,
    J1p,
    dims::NTuple{3, Int} = (4, 4, 1),
    seed::Integer = 0,
)
    all(>(0), dims) || throw(ArgumentError("all supercell dimensions must be positive"))
    exchanges_eV = (; J0, J1, J2, J3, J1p)
    all(isfinite, values(exchanges_eV)) ||
        throw(ArgumentError("all exchange constants must be finite"))

    lattice_constants = (3.85, 3.85, 20.0)
    positions = [[0.0, 0.0, 0.25], [0.0, 0.0, 0.75]]
    latvecs = lattice_vectors(lattice_constants..., 90, 90, 90)
    cryst = Crystal(latvecs, positions; types = ["Ni_top", "Ni_bottom"])
    sys = System(
        cryst,
        [1 => Moment(s = 1 / 2, g = 2), 2 => Moment(s = 1 / 2, g = 2)],
        :dipole;
        dims,
        seed,
    )
    bonds = _sunny_bilayer_bonds()
    for exchange_name in keys(bonds)
        exchange_meV = SUNNY_ENERGY_SCALE * getproperty(exchanges_eV, exchange_name)
        for (atom1, atom2, offset) in getproperty(bonds, exchange_name)
            set_exchange!(sys, exchange_meV, Bond(atom1, atom2, offset))
        end
    end
    return (;
        sys,
        cryst,
        lattice_constants,
        positions,
        dims,
        seed = Int(seed),
        spin = 0.5,
        g = 2.0,
        input_unit = "eV",
        sunny_unit = "meV",
        exchanges_eV,
        bonds,
    )
end

function _spins_matrix(sys)
    spins = zeros(Float64, length(sys.dipoles), 3)
    for (index, spin) in enumerate(sys.dipoles)
        spins[index, :] .= Tuple(spin)
    end
    return spins
end

function minimize_sunny_bilayer(; kwargs...)
    model = build_sunny_bilayer(; kwargs...)
    randomize_spins!(model.sys)
    initial_energy_eV = energy(model.sys) / SUNNY_ENERGY_SCALE
    optimization = minimize_energy!(model.sys)
    final_energy_eV = energy(model.sys) / SUNNY_ENERGY_SCALE
    spins = _spins_matrix(model.sys)
    all(isfinite, spins) || error("Sunny returned non-finite spins")
    isfinite(initial_energy_eV) || error("Sunny returned non-finite initial energy")
    isfinite(final_energy_eV) || error("Sunny returned non-finite final energy")
    phase = classify_ground_state(spins, model.dims)
    return (;
        model,
        sys = model.sys,
        initial_energy_eV,
        final_energy_eV,
        energy_per_site_eV = final_energy_eV / size(spins, 1),
        converged = hasproperty(optimization, :converged) ?
            optimization.converged : true,
        spins,
        phase,
        layer_correlation = phase.layer_correlation,
    )
end

function ground_state_structure_factor(spins::AbstractMatrix, dims::NTuple{3, Int})
    Lx, Ly, Lz = dims
    Lz == 1 || throw(ArgumentError("only one bilayer repeat along z is supported"))
    n_cells = Lx * Ly
    size(spins) == (2n_cells, 3) ||
        throw(DimensionMismatch("spins must contain two sites per in-plane cell"))
    weights = zeros(Float64, Lx, Ly)
    for mx in 0:(Lx - 1), my in 0:(Ly - 1)
        amplitude = zeros(ComplexF64, 3)
        for layer in 0:1, y in 0:(Ly - 1), x in 0:(Lx - 1)
            index = layer * n_cells + y * Lx + x + 1
            phase = exp(-2pi * im * (mx * x / Lx + my * y / Ly))
            amplitude .+= phase .* spins[index, :]
        end
        weights[mx + 1, my + 1] = sum(abs2, amplitude) / (2n_cells)
    end
    maximum_index = argmax(weights)
    mx, my = Tuple(maximum_index)
    q_rlu = [(mx - 1) / Lx, (my - 1) / Ly, 0.0]
    return (; weights, q_rlu, peak_weight = weights[maximum_index])
end

function _periodic_distance(value, target)
    delta = abs(value - target)
    return min(delta, abs(delta - 1))
end

function _matches_q(q, target; tolerance = 1e-8)
    return _periodic_distance(q[1], target[1]) <= tolerance &&
        _periodic_distance(q[2], target[2]) <= tolerance
end

function classify_ground_state(spins::AbstractMatrix, dims::NTuple{3, Int})
    structure = ground_state_structure_factor(spins, dims)
    q = structure.q_rlu
    phase = if _matches_q(q, (0.0, 0.0))
        "FM"
    elseif _matches_q(q, (0.5, 0.5))
        "Neel"
    elseif _matches_q(q, (0.5, 0.0)) || _matches_q(q, (0.0, 0.5))
        "stripe"
    elseif any(_periodic_distance(component, 0.25) <= 1e-8 ||
               _periodic_distance(component, 0.75) <= 1e-8 for component in q[1:2])
        "period-4"
    else
        "other"
    end
    n_cells = dims[1] * dims[2]
    layer_correlation = sum(
        dot(spins[index, :], spins[n_cells + index, :]) for index in 1:n_cells
    ) / (n_cells * 0.5^2)
    magnetization = vec(sum(spins; dims = 1)) / size(spins, 1)
    return (;
        phase,
        q_rlu = q,
        peak_weight = structure.peak_weight,
        structure_factor = structure.weights,
        layer_correlation,
        magnetization,
    )
end
