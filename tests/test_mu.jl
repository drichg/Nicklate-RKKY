using Test
using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "NickelateRKKY.jl"))
using .NickelateRKKY

@testset "chemical potential filling" begin
	Nk = 20
	U = 0.0
	temperature = 0.01
	n_target = 1.5

	mu = find_chemical_potential(n_target, U; Nk, temperature)
	n = filling_number(mu, U; Nk, temperature)
	n_per_orbital = filling_number_per_orbital(mu, U; Nk, temperature)

	@test isfinite(mu)
	@test isapprox(n, n_target; atol = 1e-8)
	@test length(n_per_orbital) == 4
	@test isapprox(sum(n_per_orbital), n; atol = 1e-8)
	@test all(>(0), n_per_orbital)
end
