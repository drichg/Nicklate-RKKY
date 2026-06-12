using Test

include(joinpath(
	@__DIR__,
	"..",
	"scripts",
	"plotting",
	"plot_chi_highsym.jl",
))

@testset "orbital input and projection" begin
	@test normalize_orbital("x2") == "x2"
	@test normalize_orbital("X2") == "x2"
	@test normalize_orbital(:z2) == "z2"
	@test_throws ArgumentError normalize_orbital("xy")

	chi = reshape(ComplexF64.(1:16), 4, 4)
	@test orbital_chi_value(chi, "x2") == real(chi[1, 1] + chi[1, 3])
	@test orbital_chi_value(chi, "z2") == real(chi[2, 2] + chi[2, 4])
end

@testset "simple high-symmetry path" begin
	qpoints, distances, ticks, labels = high_symmetry_path(4)

	@test isapprox(qpoints[1, :], [0.0, 0.0])
	@test isapprox(qpoints[3, :], [0.0, pi])
	@test isapprox(qpoints[5, :], [pi, pi])
	@test isapprox(qpoints[6, :], [pi / 2, pi / 2])
	@test isapprox(qpoints[end, :], [0.0, 0.0])
	@test ticks == distances[[1, 3, 5, 6, 7]]
	@test labels == ["Gamma", "Y", "M", "P", "Gamma"]
	@test_throws ArgumentError high_symmetry_path(6)
end

@testset "small cached orbital calculation and plot" begin
	tmpdir = mktempdir()
	outfile = joinpath(tmpdir, "chi_x2_test.png")

	x2_data = compute_chi_path(0.0, "x2"; Nk = 4, temperature = 0.01)
	z2_data = compute_chi_path(0.0, :z2; Nk = 4, temperature = 0.01)
	@test x2_data.orbital == "x2"
	@test z2_data.orbital == "z2"
	@test x2_data.cache.Nk == 4
	@test length(x2_data.values) == 7
	@test all(isfinite, x2_data.values)
	@test all(isfinite, z2_data.values)

	chi00 = chi_matrix_q_index(0, 0, x2_data.cache, x2_data.mu; temperature = 0.01)
	@test isapprox(x2_data.values[1], real(chi00[1, 1] + chi00[1, 3]); atol = 1e-12)
	@test isapprox(z2_data.values[1], real(chi00[2, 2] + chi00[2, 4]); atol = 1e-12)

	plt, plotted = plot_chi_path(
		0.0,
		"x2";
		Nk = 4,
		temperature = 0.01,
		savepath = outfile,
	)

	@test plt !== nothing
	@test plotted.values == x2_data.values
	@test isfile(outfile)
	@test filesize(outfile) > 0
end
