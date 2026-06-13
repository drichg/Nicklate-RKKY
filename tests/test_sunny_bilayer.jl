using Test

@testset "Sunny dependency API" begin
    using Sunny

    @test isdefined(Sunny, :System)
    @test isdefined(Sunny, :set_exchange!)
    @test isdefined(Sunny, :randomize_spins!)
    @test isdefined(Sunny, :minimize_energy!)
    @test isdefined(Sunny, :SpinWaveTheory)
    @test isdefined(Sunny, :ssf_perp)
    @test isdefined(Sunny, :intensities_bands)
end

include(joinpath(@__DIR__, "..", "src", "NickelateRKKY.jl"))
using .NickelateRKKY

@testset "vertically registered bilayer geometry" begin
    model = build_sunny_bilayer(
        J0 = 0.001,
        J1 = 0.002,
        J2 = -0.0003,
        J3 = 0.0004,
        J1p = 0.0005,
        dims = (4, 4, 1),
        seed = 11,
    )
    @test model.positions == [[0.0, 0.0, 0.25], [0.0, 0.0, 0.75]]
    @test model.lattice_constants == (3.85, 3.85, 20.0)
    @test model.dims == (4, 4, 1)
    @test length(model.sys.dipoles) == 32
    @test model.input_unit == "eV"
    @test model.sunny_unit == "meV"
    @test model.exchanges_eV.J1p == 0.0005
    @test Set(keys(model.bonds)) == Set((:J0, :J1, :J2, :J3, :J1p))
    @test model.bonds.J0 == [(1, 2, (0, 0, 0))]
    @test (1, 2, (1, 0, 0)) in model.bonds.J1p
end

@testset "deterministic Sunny minimization" begin
    kwargs = (
        J0 = 0.001,
        J1 = 0.002,
        J2 = -0.0003,
        J3 = 0.0004,
        J1p = 0.0005,
        dims = (4, 4, 1),
        seed = 20260613,
    )
    first_result = minimize_sunny_bilayer(; kwargs...)
    second_result = minimize_sunny_bilayer(; kwargs...)
    @test isfinite(first_result.initial_energy_eV)
    @test isfinite(first_result.final_energy_eV)
    @test first_result.final_energy_eV <= first_result.initial_energy_eV + 1e-10
    @test size(first_result.spins) == (32, 3)
    @test all(isapprox.(sqrt.(sum(abs2, first_result.spins; dims = 2)), 0.5; atol = 1e-10))
    @test first_result.final_energy_eV ≈ second_result.final_energy_eV atol = 1e-12
    @test first_result.spins ≈ second_result.spins atol = 1e-10
end

function synthetic_spins(Lx, Ly, pattern)
    spins = zeros(Float64, Lx * Ly * 2, 3)
    index = 1
    for layer in 1:2, y in 0:(Ly - 1), x in 0:(Lx - 1)
        sign = pattern(x, y, layer)
        spins[index, 3] = 0.5 * sign
        index += 1
    end
    return spins
end

@testset "ground-state phase classification" begin
    @test classify_ground_state(
        synthetic_spins(4, 4, (x, y, layer) -> 1),
        (4, 4, 1),
    ).phase == "FM"
    @test classify_ground_state(
        synthetic_spins(4, 4, (x, y, layer) -> (-1)^(x + y)),
        (4, 4, 1),
    ).phase == "Neel"
    @test classify_ground_state(
        synthetic_spins(4, 4, (x, y, layer) -> (-1)^x),
        (4, 4, 1),
    ).phase == "stripe"
    @test classify_ground_state(
        synthetic_spins(4, 4, (x, y, layer) -> x < 2 ? 1 : -1),
        (4, 4, 1),
    ).phase == "period-4"
end

@testset "continuous classical exchange minimum" begin
    ferromagnet = classical_bilayer_minimum(
        (; J0 = 0.0, J1 = -1.0, J2 = 0.0, J3 = 0.0, J1p = 0.0);
        grid_size = 101,
    )
    neel = classical_bilayer_minimum(
        (; J0 = 0.0, J1 = 1.0, J2 = 0.0, J3 = 0.0, J1p = 0.0);
        grid_size = 101,
    )
    period4 = classical_bilayer_minimum(
        (; J0 = 1.0, J1 = 0.0, J2 = 0.0, J3 = 1.0, J1p = 0.0);
        grid_size = 101,
    )
    @test ferromagnet.q_rlu[1:2] ≈ [0.0, 0.0] atol = 1e-12
    @test neel.q_rlu[1:2] ≈ [0.5, 0.5] atol = 1e-12
    @test period4.q_rlu[1:2] ≈ [0.25, 0.25] atol = 0.005
    @test period4.layer_parity == -1
end

@testset "Sunny path and instability classification" begin
    model = build_sunny_bilayer(
        J0 = 0.001,
        J1 = 0.002,
        J2 = 0.0,
        J3 = 0.0,
        J1p = 0.0,
        dims = (1, 1, 1),
        seed = 1,
    )
    path = sunny_bilayer_q_path(model.cryst; n_points = 301)
    @test size(path.q_points) == (301, 3)
    @test path.q_points[1, :] == [0.0, 0.0, 0.0]
    @test path.q_points[end, :] == [0.0, 0.0, 0.0]
    @test path.labels == ["Gamma", "X", "M", "Gamma"]
    @test issorted(path.distances)

    first_error = ErrorException("Instability at wavevector q = [0.5, 0.0, 0.0]")
    classified = classify_sunny_exception(first_error; stage = "spectrum")
    @test classified.status == "unstable_spin_wave"
    @test classified.stage == "spectrum"
    @test classified.q == [0.5, 0.0, 0.0]

    second_error = ErrorException("Not an energy-minimum; wavevector q = [0, 0, 0] unstable.")
    @test classify_sunny_exception(second_error).status == "unstable_spin_wave"
    @test classify_sunny_exception(ErrorException("file missing")).status == "software_error"
end

@testset "stable retry selection" begin
    attempts = [
        (; status = "unstable_spin_wave", final_energy_eV = -1.0),
        (; status = "ok", final_energy_eV = -0.8),
        (; status = "ok", final_energy_eV = -1.2),
    ]
    @test select_sunny_attempt(attempts) == 3
    @test select_sunny_attempt(attempts[1:1]) == 1
end

@testset "Sunny spectrum and weighted plot" begin
    minimized = minimize_sunny_bilayer(
        J0 = 0.001,
        J1 = -0.002,
        J2 = 0.0,
        J3 = 0.0,
        J1p = 0.0,
        dims = (2, 2, 1),
        seed = 4,
    )
    spectrum = compute_sunny_bilayer_spectrum(minimized.sys; n_q = 31)
    @test spectrum.status == "ok"
    @test size(spectrum.energies_eV) == size(spectrum.weights)
    @test size(spectrum.energies_eV, 2) == 31
    @test all(isfinite, spectrum.energies_eV)
    @test all(isfinite, spectrum.weights)
    @test all(spectrum.weights .>= -1e-12)
    @test spectrum.measure == "ssf_perp"
    @test spectrum.energy_unit == "eV"

    stable = classify_spin_wave([0.0 -1e-10; 0.2 0.3]; negative_tolerance = 1e-8)
    unstable = classify_spin_wave([0.0 -0.01; 0.2 0.3]; negative_tolerance = 1e-8)
    @test stable.status == "ok"
    @test unstable.status == "unstable_spin_wave"
    @test unstable.negative_count == 1
    @test unstable.minimum_frequency_eV == -0.01

    mktempdir() do tmpdir
        filepath = joinpath(tmpdir, "spectrum.png")
        plot_sunny_bilayer_spectrum(spectrum; U = 0.0, filepath)
        @test isfile(filepath)
        @test filesize(filepath) > 0
    end
end

@testset "Sunny retry protocol" begin
    calls = Tuple[]
    attempt_calculator = function(U, exchanges; dims, seed, kwargs...)
        push!(calls, (dims, seed))
        status = dims == (6, 6, 1) ? "ok" : "unstable_spin_wave"
        return (;
            status,
            final_energy_eV = dims == (6, 6, 1) ? -seed / 1e9 : -1.0,
            dims,
            seed,
        )
    end
    retried = run_sunny_point_with_retries(
        2.0,
        (; J0 = 0.0, J1 = 0.0, J2 = 0.0, J3 = 0.0, J1p = 0.0);
        initial_seed = 100,
        attempt_calculator,
    )
    @test length(retried.attempts) == 7
    @test retried.selected_attempt in 6:7
    @test retried.status == "ok"
    @test all(call[1] != (8, 8, 1) for call in calls)
end

@testset "atomic Sunny persistence and resume" begin
    rkky = (;
        metadata = (; schema_version = 1, exchange_unit = "eV"),
        completed_U_values = [0.0, 0.1],
        J0 = [0.1, 0.2],
        J1 = [0.2, 0.3],
        J2 = [0.3, 0.4],
        J3 = [0.4, 0.5],
        J1p = [0.5, 0.6],
    )
    mktempdir() do tmpdir
        filepath = joinpath(tmpdir, "sunny.jld2")
        image_dir = joinpath(tmpdir, "images")
        calls = Float64[]
        fail_once = Ref(true)
        calculator = function(U, exchanges; initial_seed, kwargs...)
            push!(calls, U)
            if U == 0.1 && fail_once[]
                fail_once[] = false
                error("simulated interruption")
            end
            return (;
                U,
                status = "ok",
                selected_attempt = 1,
                attempts = [(;
                    status = "ok",
                    final_energy_eV = -U,
                    energy_per_site_eV = -U,
                    phase = (; phase = "FM", q_rlu = [0.0, 0.0, 0.0]),
                    spectrum = (;
                        status = "ok",
                        energies_eV = [0.0 0.1],
                        weights = [1.0 1.0],
                        distances = [0.0, 1.0],
                        tick_distances = [0.0, 1.0],
                        labels = ["Gamma", "Gamma"],
                    ),
                )],
            )
        end
        plot_writer = function(point, image_path)
            write(image_path, "plot")
        end

        @test_throws ErrorException scan_sunny_bilayer(
            rkky;
            output_path = filepath,
            image_dir,
            point_calculator = calculator,
            plot_writer,
        )
        partial = load_sunny_result(filepath)
        @test partial.completed_U_values == [0.0]
        resumed = scan_sunny_bilayer(
            rkky;
            output_path = filepath,
            image_dir,
            point_calculator = calculator,
            plot_writer,
        )
        @test resumed.completed_U_values == [0.0, 0.1]
        @test calls == [0.0, 0.1, 0.1]
        @test !isfile(filepath * ".tmp")
    end
end

@testset "default Sunny script configuration" begin
    script = joinpath(
        @__DIR__,
        "..",
        "scripts",
        "sunny",
        "run_bilayer_groundstate_spectrum.jl",
    )
    @test isfile(script)
    include(script)
    config = default_sunny_bilayer_config()
    @test config.dims == (4, 4, 1)
    @test config.base_seed == 20260612
    @test config.n_q == 301
    @test endswith(config.input_path, "rkky_x2_bilayer_U0-6_Nk100.jld2")
    @test endswith(config.output_path, "sunny_bilayer_U0-6.jld2")
    @test sunny_runtime_version() == "0.8.0"
end
