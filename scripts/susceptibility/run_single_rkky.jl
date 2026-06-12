using LinearAlgebra
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))
using Plots

include(joinpath(@__DIR__, "..", "..", "src", "NickelateRKKY.jl"))
using .NickelateRKKY

if abspath(PROGRAM_FILE) == @__FILE__
let
	Nk = 50
	Nq = Nk
	U = 0.0
	n_target = 1.5
	temperature = 0.001

	mu = find_chemical_potential(n_target, U; Nk = Nk, temperature = temperature)
	println("mu = ", mu)
	println("n = ", filling_number(mu, U; Nk = Nk, temperature = temperature))

	cache = diagonalized_kmesh(Nk, U = U)



	J1 = 0.0 + 0.0im
	J2 = 0.0 + 0.0im
	J3 = 0.0 + 0.0im
	J0 = 0.0 + 0.0im
	for iqx in 0:(Nk-1), iqy in 0:(Nk-1)
		qx = 2pi * iqx / Nk
		qy = 2pi * iqy / Nk
		chi = chi_matrix_q_index(iqx, iqy, cache, mu; temperature = temperature)
		# AA layer, dx2-y2 channel
		# chiq = chi[2, 4]
		# chiq=sum(chi)


		chiq = chi[2, 2]
		J1 += exp(im * qx) * chiq          # R = (1,0)
		J2 += exp(im * (qx + qy)) * chiq   # R = (1,1)
		J3 += exp(im * 2qx) * chiq         # R = (2,0)
		J0 += chiq
	end

	J1 = -real(J1) / Nk^2
	J2 = -real(J2) / Nk^2
	J3 = -real(J3) / Nk^2
	J0 = -real(J0) / Nk^2
	println("J1 = ", J1)
	println("J2 = ", J2)
	println("J3 = ", J3)
	println("J0 = ", J0)
	println("|J1|:|J2|:|J3|:|J0| = ", abs(J1), " : ", abs(J2), " : ", abs(J3), " : ", abs(J0))
	println("|J3/J1| = ", abs(J3 / J1))
end
end
