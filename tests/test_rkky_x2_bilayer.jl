using JLD2
using Test

include(joinpath(@__DIR__, "..", "src", "NickelateRKKY.jl"))
using .NickelateRKKY

@testset "x2 layer projection" begin
    chi = reshape(ComplexF64.(1:16), 4, 4)
    channels = x2_layer_channels(chi)

    @test X1_INDEX == 1
    @test X2_INDEX == 3
    @test channels.intra == (chi[1, 1] + chi[3, 3]) / 2
    @test channels.inter == (chi[1, 3] + chi[3, 1]) / 2
    @test channels.intra_layer_residual == abs(chi[1, 1] - chi[3, 3])
    @test channels.inter_layer_residual == abs(chi[1, 3] - chi[3, 1])
    @test_throws DimensionMismatch x2_layer_channels(zeros(3, 3))
end

@testset "five x2 bilayer exchanges" begin
    Nk = 4
    U = 0.0
    temperature = 0.01
    n_target = 1.5
    JK = 2.0
    mu = find_chemical_potential(n_target, U; Nk, temperature)
    cache = diagonalized_kmesh(Nk, U)
    sums = zeros(ComplexF64, 5)

    for iqx in 0:(Nk - 1), iqy in 0:(Nk - 1)
        qx = 2pi * iqx / Nk
        qy = 2pi * iqy / Nk
        channels = x2_layer_channels(
            chi_matrix_q_index(iqx, iqy, cache, mu; temperature),
        )
        sums[1] += channels.inter
        sums[2] += exp(im * qx) * channels.intra
        sums[3] += exp(im * (qx + qy)) * channels.intra
        sums[4] += exp(im * 2qx) * channels.intra
        sums[5] += exp(im * qx) * channels.inter
    end

    bare = -real.(sums) ./ Nk^2
    actual = compute_x2_bilayer_rkky(
        U;
        Nk,
        Nq = Nk,
        n_target,
        temperature,
        JK,
    )

    @test isapprox(
        [actual.bare_J0, actual.bare_J1, actual.bare_J2, actual.bare_J3, actual.bare_J1p],
        bare;
        atol = 1e-12,
    )
    @test isapprox(
        [actual.J0, actual.J1, actual.J2, actual.J3, actual.J1p],
        JK^2 .* bare;
        atol = 1e-12,
    )
    @test isapprox(actual.filling, n_target; atol = 1e-8)
    @test actual.max_fourier_imaginary_residual < 1e-10
    @test_throws ArgumentError compute_x2_bilayer_rkky(U; Nk = 4, Nq = 5)
end

@testset "RKKY metadata and resume validation" begin
    U_values = collect(0.0:0.1:0.2)
    metadata = x2_rkky_metadata(
        U_values;
        Nk = 4,
        Nq = 4,
        n_target = 1.5,
        temperature = 0.01,
        JK = 1.0,
    )

    @test metadata.schema_version == 1
    @test metadata.basis == "x1,z1,x2,z2"
    @test metadata.exchange_unit == "eV"
    @test metadata.susceptibility_unit == "eV^-1"
    @test metadata.U_values == U_values
    @test validate_x2_rkky_resume(metadata, metadata) === nothing

    changed = merge(metadata, (; Nk = 6))
    @test_throws ArgumentError validate_x2_rkky_resume(metadata, changed)
end

@testset "atomic RKKY scan persistence" begin
    mktempdir() do tmpdir
        filepath = joinpath(tmpdir, "rkky.jld2")
        calls = Float64[]
        fail_once = Ref(true)
        calculator = function(U; Nk, Nq, n_target, temperature, JK)
            push!(calls, U)
            if U == 0.1 && fail_once[]
                fail_once[] = false
                error("simulated interruption")
            end
            return (;
                U,
                mu = U + 1,
                filling = n_target,
                filling_error = 0.0,
                bare_J0 = U + 0.1,
                bare_J1 = U + 0.2,
                bare_J2 = U + 0.3,
                bare_J3 = U + 0.4,
                bare_J1p = U + 0.5,
                J0 = JK^2 * (U + 0.1),
                J1 = JK^2 * (U + 0.2),
                J2 = JK^2 * (U + 0.3),
                J3 = JK^2 * (U + 0.4),
                J1p = JK^2 * (U + 0.5),
                max_intra_layer_residual = 0.0,
                max_inter_layer_residual = 0.0,
                max_fourier_imaginary_residual = 0.0,
            )
        end

        @test_throws ErrorException scan_x2_bilayer_rkky(
            [0.0, 0.1];
            Nk = 4,
            Nq = 4,
            filepath,
            point_calculator = calculator,
        )
        @test isfile(filepath)
        first_result = load_x2_rkky_result(filepath)
        @test first_result.completed_U_values == [0.0]
        @test calls == [0.0, 0.1]
        @test !isfile(filepath * ".tmp")

        resumed = scan_x2_bilayer_rkky(
            [0.0, 0.1];
            Nk = 4,
            Nq = 4,
            filepath,
            point_calculator = calculator,
        )
        @test resumed.completed_U_values == [0.0, 0.1]
        @test calls == [0.0, 0.1, 0.1]
        @test load_x2_rkky_result(filepath) == resumed
    end
end

@testset "default RKKY script configuration" begin
    script = joinpath(
        @__DIR__,
        "..",
        "scripts",
        "susceptibility",
        "run_x2_bilayer_rkky.jl",
    )
    @test isfile(script)
    include(script)
    config = default_x2_bilayer_rkky_config()
    @test config.U_values == collect(0.0:0.1:6.0)
    @test config.Nk == 100
    @test config.Nq == 100
    @test config.n_target == 1.5
    @test config.temperature == 0.01
    @test config.JK == 1.0
    @test endswith(config.filepath, "rkky_x2_bilayer_U0-6_Nk100.jld2")
end

@testset "RKKY convergence script configuration" begin
    script = joinpath(
        @__DIR__,
        "..",
        "scripts",
        "susceptibility",
        "check_x2_rkky_convergence.jl",
    )
    @test isfile(script)
    include(script)
    config = default_x2_convergence_config()
    @test config.U_values == [0.0, 2.9, 3.0, 3.6, 6.0]
    @test config.base_grid == 80
    @test config.reference_grid == 100
    @test config.diagnostic_grid == 120
    @test config.relative_tolerance == 0.05
end
