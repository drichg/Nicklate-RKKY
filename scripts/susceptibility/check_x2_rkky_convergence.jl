if !isdefined(Main, :NickelateRKKY)
    include(joinpath(@__DIR__, "..", "..", "src", "NickelateRKKY.jl"))
end
using .NickelateRKKY
using JLD2

function default_x2_convergence_config()
    return (;
        U_values = [0.0, 2.9, 3.0, 3.6, 6.0],
        base_grid = 80,
        reference_grid = 100,
        diagnostic_grid = 120,
        relative_tolerance = 0.05,
        absolute_tolerance_eV = 1e-6,
        input_path = joinpath(
            NickelateRKKY.PROJECT_ROOT,
            "data",
            "results",
            "rkky_x2_bilayer_U0-6_Nk100.jld2",
        ),
        output_path = joinpath(
            NickelateRKKY.PROJECT_ROOT,
            "data",
            "results",
            "rkky_x2_bilayer_convergence.jld2",
        ),
    )
end

function _exchange_vector(point)
    return [point.J0, point.J1, point.J2, point.J3, point.J1p]
end

function run_x2_rkky_convergence()
    config = default_x2_convergence_config()
    reference = load_x2_rkky_result(config.input_path)
    names = ["J0", "J1", "J2", "J3", "J1p"]
    rows = Any[]
    for U in config.U_values
        index = findfirst(==(U), reference.completed_U_values)
        index === nothing && error("U=$U is missing from the reference scan")
        reference_values = [
            reference.J0[index],
            reference.J1[index],
            reference.J2[index],
            reference.J3[index],
            reference.J1p[index],
        ]
        base = compute_x2_bilayer_rkky(
            U;
            Nk = config.base_grid,
            Nq = config.base_grid,
            n_target = 1.5,
            temperature = 0.01,
            JK = 1.0,
        )
        base_values = _exchange_vector(base)
        tolerances = max.(
            config.relative_tolerance .* abs.(reference_values),
            config.absolute_tolerance_eV,
        )
        differences = abs.(base_values .- reference_values)
        needs_diagnostic = any(differences .> tolerances)
        diagnostic = needs_diagnostic ? compute_x2_bilayer_rkky(
            U;
            Nk = config.diagnostic_grid,
            Nq = config.diagnostic_grid,
            n_target = 1.5,
            temperature = 0.01,
            JK = 1.0,
        ) : nothing
        push!(rows, (;
            U,
            names,
            grid80 = base_values,
            grid100 = reference_values,
            absolute_difference_80_100 = differences,
            tolerances,
            needs_diagnostic,
            grid120 = diagnostic === nothing ? nothing : _exchange_vector(diagnostic),
        ))
        println("U=$U, diagnostic Nk=120: $needs_diagnostic")
    end
    result = (; config, rows)
    mkpath(dirname(config.output_path))
    jldsave(config.output_path; result)
    println("Saved convergence result: ", config.output_path)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_x2_rkky_convergence()
end
