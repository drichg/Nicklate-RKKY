if !isdefined(Main, :NickelateRKKY)
    include(joinpath(@__DIR__, "..", "..", "src", "NickelateRKKY.jl"))
end
using .NickelateRKKY

function default_x2_bilayer_rkky_config()
    return (;
        U_values = collect(0.0:0.1:6.0),
        Nk = 100,
        Nq = 100,
        n_target = 1.5,
        temperature = 0.01,
        JK = 1.0,
        filepath = joinpath(
            NickelateRKKY.PROJECT_ROOT,
            "data",
            "results",
            "rkky_x2_bilayer_U0-6_Nk100.jld2",
        ),
    )
end

function run_default_x2_bilayer_rkky()
    config = default_x2_bilayer_rkky_config()
    println("Starting bilayer x2 RKKY scan with configuration:")
    println(config)
    result = scan_x2_bilayer_rkky(config.U_values; (
        field => getproperty(config, field) for
        field in (:Nk, :Nq, :n_target, :temperature, :JK, :filepath)
    )...)
    println("Completed $(length(result.completed_U_values))/$(length(config.U_values)) U points")
    println("Saved result: $(config.filepath)")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_default_x2_bilayer_rkky()
end
