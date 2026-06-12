using LinearAlgebra
using DelimitedFiles

include(joinpath(@__DIR__, "..", "..", "src", "NickelateRKKY.jl"))
using .NickelateRKKY

function rkky_all_interactions(U_values;
	Nk::Integer = 50,
	n_target = 1.5,
	temperature = 0.001,
)
	Nq = Nk
	channel_labels = ["orbital1", "orbital2", "orbital_sum"]
	n_channels = length(channel_labels)
	n_U = length(U_values)

	J0 = zeros(Float64, n_U, n_channels)  # interlayer same-orbital: chi13, chi24
	J1 = zeros(Float64, n_U, n_channels)  # intralayer R = (1,0): chi11, chi22
	J2 = zeros(Float64, n_U, n_channels)  # intralayer R = (1,1): chi11, chi22
	J3 = zeros(Float64, n_U, n_channels)  # intralayer R = (2,0): chi11, chi22
	Js = [J0, J1, J2, J3]

	mu_values = zeros(Float64, n_U)
	filling_values = zeros(Float64, n_U)

	for (iU, U) in pairs(U_values)
		mu = find_chemical_potential(n_target, U; Nk = Nk, temperature = temperature)
		filling = filling_number(mu, U; Nk = Nk, temperature = temperature)
		mu_values[iU] = mu
		filling_values[iU] = filling
		println("U = ", U, ", mu = ", mu, ", n = ", filling)

		cache = diagonalized_kmesh(Nk, U)
		Jsum = zeros(ComplexF64, 4, n_channels)

		for iqx in 0:(Nq-1), iqy in 0:(Nq-1)
			qx = 2pi * iqx / Nq
			qy = 2pi * iqy / Nq
			chi = chi_matrix_q_index(iqx, iqy, cache, mu; temperature = temperature)

			interlayer_channels = (chi[1, 3], chi[2, 4], chi[1, 3] + chi[2, 4])
			intralayer_channels = (chi[1, 1], chi[2, 2], chi[1, 1] + chi[2, 2])

			for channel in 1:n_channels
				Jsum[1, channel] += interlayer_channels[channel]
				Jsum[2, channel] += exp(im * qx) * intralayer_channels[channel]
				Jsum[3, channel] += exp(im * (qx + qy)) * intralayer_channels[channel]
				Jsum[4, channel] += exp(im * 2qx) * intralayer_channels[channel]
			end
		end

		for j in 1:4
			Js[j][iU, :] .= -real.(Jsum[j, :]) ./ Nq^2
		end

		println("J0 interlayer = ", J0[iU, :])
		println("J1 intralayer R=(1,0) = ", J1[iU, :])
		println("J2 intralayer R=(1,1) = ", J2[iU, :])
		println("J3 intralayer R=(2,0) = ", J3[iU, :])
	end

	return (;
		U_values,
		mu_values,
		filling_values,
		channel_labels,
		Js,
		J0,
		J1,
		J2,
		J3,
		J_descriptions = [
			"J0: interlayer same-orbital, channels chi13, chi24, chi13+chi24",
			"J1: intralayer R=(1,0), channels chi11, chi22, chi11+chi22",
			"J2: intralayer R=(1,1), channels chi11, chi22, chi11+chi22",
			"J3: intralayer R=(2,0), channels chi11, chi22, chi11+chi22",
		],
	)
end

if abspath(PROGRAM_FILE) == @__FILE__
	U = range(0.0, stop = 8.0, step = 0.1)
	rslt = rkky_all_interactions(U)
end
