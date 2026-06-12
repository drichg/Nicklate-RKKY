using JLD2
using Test

run2_path = joinpath(
	@__DIR__,
	"..",
	"scripts",
	"susceptibility",
	"run_total_rkky.jl",
)

@testset "run2 total RKKY" begin
	@test isfile(run2_path)

	if isfile(run2_path)
		include(run2_path)

		Nk = 4
		temperature = 0.01
		U = 0.0
		rslt = rkky_total_interactions([U]; Nk = Nk, temperature = temperature)

		@test rslt.channel_labels == ["total"]
		@test size(rslt.J0) == (1,)
		@test size(rslt.J1) == (1,)
		@test size(rslt.J2) == (1,)
		@test size(rslt.J3) == (1,)
		@test isapprox(rslt.filling_values[1], 1.5; atol = 1e-8)

		cache = diagonalized_kmesh(Nk, U)
		Jsum = zeros(ComplexF64, 4)
		for iqx in 0:(Nk-1), iqy in 0:(Nk-1)
			qx = 2pi * iqx / Nk
			qy = 2pi * iqy / Nk
			chi = chi_matrix_q_index(iqx, iqy, cache, rslt.mu_values[1]; temperature)

			total_interlayer = sum(chi[1:2, 3:4])
			total_intralayer = sum(chi[1:2, 1:2])

			Jsum[1] += total_interlayer
			Jsum[2] += exp(im * qx) * total_intralayer
			Jsum[3] += exp(im * (qx + qy)) * total_intralayer
			Jsum[4] += exp(im * 2qx) * total_intralayer
		end

		expected = -real.(Jsum) ./ Nk^2
		@test isapprox(rslt.J0[1], expected[1]; atol = 1e-12)
		@test isapprox(rslt.J1[1], expected[2]; atol = 1e-12)
		@test isapprox(rslt.J2[1], expected[3]; atol = 1e-12)
		@test isapprox(rslt.J3[1], expected[4]; atol = 1e-12)

		savepath = tempname() * ".jld2"
		saved = rkky_total_interactions([U]; Nk = Nk, temperature = temperature, savepath = savepath)
		@test isfile(savepath)
		loaded = load(savepath, "rslt")
		@test loaded.U_values == saved.U_values
		@test loaded.J0 == saved.J0
		@test loaded.J1 == saved.J1
		@test loaded.J2 == saved.J2
		@test loaded.J3 == saved.J3
	end
end
