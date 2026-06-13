if !isdefined(Main, :NickelateRKKY)
    include(joinpath(@__DIR__, "..", "..", "src", "NickelateRKKY.jl"))
end
using .NickelateRKKY

function default_sunny_bilayer_config()
    return (;
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
            "sunny_bilayer_U0-6.jld2",
        ),
        image_dir = joinpath(
            NickelateRKKY.PROJECT_ROOT,
            "results",
            "excitations",
        ),
        dims = (4, 4, 1),
        base_seed = 20260612,
        n_q = 301,
        negative_tolerance = 1e-8,
    )
end

function run_default_sunny_bilayer()
    config = default_sunny_bilayer_config()
    rkky_result = load_x2_rkky_result(config.input_path)
    rkky_result.completed_U_values == collect(0.0:0.1:6.0) ||
        error("default Sunny run requires the complete U=0:0.1:6 RKKY scan")
    println("Sunny version: ", Base.pkgversion(Sunny))
    println("Configuration: ", config)
    result = scan_sunny_bilayer(
        rkky_result;
        output_path = config.output_path,
        image_dir = config.image_dir,
        base_seed = config.base_seed,
        n_q = config.n_q,
        negative_tolerance = config.negative_tolerance,
    )
    counts = Dict{String, Int}()
    for point in result.points
        counts[point.status] = get(counts, point.status, 0) + 1
    end
    println("Completed $(length(result.completed_U_values)) points")
    println("Status counts: ", counts)
    println("Saved result: ", config.output_path)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_default_sunny_bilayer()
end
