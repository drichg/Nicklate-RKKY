using LinearAlgebra
using JLD2

include(joinpath(@__DIR__, "..", "..", "src", "NickelateRKKY.jl"))
using .NickelateRKKY

function rkky_total_interactions(U_values;
	Nk::Integer = 50,
	n_target = 1.5,
	temperature = 0.001,
	savepath::Union{Nothing, AbstractString} = nothing,
)
	Nq = Nk
	n_U = length(U_values)

	J0 = zeros(Float64, n_U)  # interlayer R=(0,0), sum chi[1:2, 3:4]
	J1 = zeros(Float64, n_U)  # intralayer R=(1,0), sum chi[1:2, 1:2]
	J2 = zeros(Float64, n_U)  # intralayer R=(1,1), sum chi[1:2, 1:2]
	J3 = zeros(Float64, n_U)  # intralayer R=(2,0), sum chi[1:2, 1:2]
	Js = [J0, J1, J2, J3]

	mu_values = zeros(Float64, n_U)
	filling_values = zeros(Float64, n_U)

	function result_prefix(n_done)
		return (;
			U_values = U_values[1:n_done],
			mu_values = mu_values[1:n_done],
			filling_values = filling_values[1:n_done],
			channel_labels = ["total"],
			Js = [J0[1:n_done], J1[1:n_done], J2[1:n_done], J3[1:n_done]],
			J0 = J0[1:n_done],
			J1 = J1[1:n_done],
			J2 = J2[1:n_done],
			J3 = J3[1:n_done],
			J_descriptions = [
				"J0: total interlayer R=(0,0), sum chi[1:2, 3:4]",
				"J1: total intralayer R=(1,0), sum chi[1:2, 1:2]",
				"J2: total intralayer R=(1,1), sum chi[1:2, 1:2]",
				"J3: total intralayer R=(2,0), sum chi[1:2, 1:2]",
			],
		)
	end

	for (iU, U) in pairs(U_values)
		mu = find_chemical_potential(n_target, U; Nk = Nk, temperature = temperature)
		filling = filling_number(mu, U; Nk = Nk, temperature = temperature)
		mu_values[iU] = mu
		filling_values[iU] = filling
		println("U = ", U, ", mu = ", mu, ", n = ", filling)

		cache = diagonalized_kmesh(Nk, U)
		Jsum = zeros(ComplexF64, 4)

		for iqx in 0:(Nq-1), iqy in 0:(Nq-1)
			qx = 2pi * iqx / Nq
			qy = 2pi * iqy / Nq
			chi = chi_matrix_q_index(iqx, iqy, cache, mu; temperature = temperature)

			total_interlayer = sum(chi[1:2, 3:4])
			total_intralayer = sum(chi[1:2, 1:2])

			Jsum[1] += total_interlayer
			Jsum[2] += exp(im * qx) * total_intralayer
			Jsum[3] += exp(im * (qx + qy)) * total_intralayer
			Jsum[4] += exp(im * 2qx) * total_intralayer
		end

		for j in 1:4
			Js[j][iU] = -real(Jsum[j]) / Nq^2
		end

		println("J0 total interlayer R=(0,0) = ", J0[iU])
		println("J1 total intralayer R=(1,0) = ", J1[iU])
		println("J2 total intralayer R=(1,1) = ", J2[iU])
		println("J3 total intralayer R=(2,0) = ", J3[iU])

		if savepath !== nothing
			jldsave(savepath; rslt = result_prefix(iU))
			println("saved partial result to ", savepath, " (", iU, "/", n_U, " U points)")
		end
	end

	return result_prefix(n_U)
end

if abspath(PROGRAM_FILE) == @__FILE__
	U = range(0.0, stop = 8, step = 0.1)
	savepath = joinpath(NickelateRKKY.PROJECT_ROOT, "data", "results", "U=0-8-total.jld2")
	mkpath(dirname(savepath))
	rslt = rkky_total_interactions(U; savepath)
end
