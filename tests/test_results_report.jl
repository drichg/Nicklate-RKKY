using Test

include(joinpath(
    @__DIR__,
    "..",
    "scripts",
    "report",
    "generate_rkky_sunny_report.jl",
))

@testset "results report configuration" begin
    config = default_results_report_config()
    @test endswith(config.rkky_path, "rkky_x2_bilayer_U0-6_Nk100.jld2")
    @test endswith(config.sunny_path, "sunny_bilayer_U0-6.jld2")
    @test endswith(config.convergence_path, "rkky_x2_bilayer_convergence.jld2")
    @test endswith(config.tex_path, "rkky_sunny_results.tex")
    @test endswith(config.figure_dir, joinpath("results", "summary"))
end
