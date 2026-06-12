# function tight_binding(kx, ky, U, mu)
# 	U = 0
# 	mu = 0
# 	t1x = -0.483
# 	t1z = -0.11
# 	t2x = 0.069
# 	t2z = -0.017
# 	t3xz = 0.239
# 	t0x = 0.005
# 	t0z = -0.635
# 	t4xz = -0.034
# 	ϵx = 0.776
# 	ϵz = 0.409

# 	Tx = 2 * t1x * (cos(kx) + cos(ky)) + 4 * t2x * cos(kx) * cos(ky) + ϵx
# 	Tz = 2 * t1z * (cos(kx) + cos(ky)) + 4 * t2z * cos(kx) * cos(ky) + ϵz
# 	Vk = 2 * t3xz * (cos(kx) - cos(ky))
# 	Vk2 = 2 * t4xz * (cos(kx) - cos(ky))

# 	HA = [Tx Vk; Vk Tz]
# 	HAB = [t0x Vk2; Vk2 t0z]

# 	return [HA HAB; HAB' HA]
# end



function fermi(energy, mu, temperature)
	temperature > 0 || throw(ArgumentError("temperature must be positive"))

	x = (energy - mu) / temperature
	if x > 40
		return 0.0
	elseif x < -40
		return 1.0
	end

	return 1 / (exp(x) + 1)
end

function fermi_derivative(energy, mu, temperature)
	temperature > 0 || throw(ArgumentError("temperature must be positive"))

	x = (energy - mu) / temperature
	if abs(x) > 80
		return 0.0
	end

	return -1 / (4 * temperature * cosh(x / 2)^2)
end

function susceptibility_kernel(e1, e2, mu, temperature; tol = 1e-8)
	delta = e1 - e2

	if abs(delta) < tol
		return -fermi_derivative((e1 + e2) / 2, mu, temperature)
	end

	return -(fermi(e1, mu, temperature) - fermi(e2, mu, temperature)) / delta
end

function band_energies_kmesh(Nk::Integer, U)
	Nk > 0 || throw(ArgumentError("Nk must be positive"))

	energies = zeros(Float64, Nk, Nk, 4)
	for ix in 1:Nk, iy in 1:Nk
		kx = 2pi * (ix - 1) / Nk
		ky = 2pi * (iy - 1) / Nk
		energies[ix, iy, :] = eigvals(Hermitian(tight_binding(kx, ky, U, 0.0)))
	end

	return energies
end

function filling_number(mu, U; Nk::Integer = 50, temperature = 0.001)
	temperature > 0 || throw(ArgumentError("temperature must be positive"))


	energies = band_energies_kmesh(Nk, U)
	return sum(fermi(energy, mu, temperature) for energy in energies) / size(energies, 1)^2
end
function filling_number_per_orbital(mu, U; Nk::Integer = 50, temperature = 0.001)
	temperature > 0 || throw(ArgumentError("temperature must be positive"))

	occupation = zeros(Float64, 4)
	for ix in 1:Nk, iy in 1:Nk
		kx = 2pi * (ix - 1) / Nk
		ky = 2pi * (iy - 1) / Nk
		eigen_k = eigen(Hermitian(tight_binding(kx, ky, U, 0.0)))

		for band in 1:4, alpha in 1:4
			occupation[alpha] += abs2(eigen_k.vectors[alpha, band]) *
								 fermi(eigen_k.values[band], mu, temperature)
		end
	end

	return occupation / Nk^2
end

function find_chemical_potential(
	n_target,
	U;
	Nk::Integer = 50,
	temperature = 0.001,
	tol = 1e-10,
	maxiter::Integer = 200,
	bracket_width = 1.0)
	Nk > 0 || throw(ArgumentError("Nk must be positive"))
	temperature > 0 || throw(ArgumentError("temperature must be positive"))


	energies = band_energies_kmesh(Nk, U)
	n_of_mu(mu) = sum(fermi(energy, mu, temperature) for energy in energies) / Nk^2
	mu0 = chemical_potential(U)
	# mu0 = 0
	width = bracket_width
	lo = mu0 - width
	hi = mu0 + width

	for _ in 1:80
		n_lo = n_of_mu(lo)
		n_hi = n_of_mu(hi)
		if n_lo <= n_target && n_hi >= n_target
			break
		end
		width *= 2
		lo = mu0 - width
		hi = mu0 + width
	end

	n_of_mu(lo) <= n_target ||
		throw(ArgumentError("failed to bracket chemical potential below n_target"))
	n_of_mu(hi) >= n_target ||
		throw(ArgumentError("failed to bracket chemical potential above n_target"))

	for _ in 1:maxiter
		mid = (lo + hi) / 2
		n_mid = n_of_mu(mid)

		if abs(n_mid - n_target) < tol || abs(hi - lo) < tol
			return mid
		elseif n_mid < n_target
			lo = mid
		else
			hi = mid
		end
	end

	return (lo + hi) / 2
end

function chi_matrix_q(qx, qy, U, mu; Nk::Integer = 50, temperature = 0.001, tol = 1e-8)
	Nk > 0 || throw(ArgumentError("Nk must be positive"))

	chi = zeros(ComplexF64, 4, 4)

	for ix in 0:(Nk-1), iy in 0:(Nk-1)
		kx = 2pi * ix / Nk
		ky = 2pi * iy / Nk

		eigen_k = eigen(Hermitian(tight_binding(kx, ky, U, 0.0)))
		eigen_kq = eigen(Hermitian(tight_binding(kx + qx, ky + qy, U, 0.0)))

		for n in 1:4, m in 1:4
			kernel = susceptibility_kernel(
				eigen_k.values[n],
				eigen_kq.values[m],
				mu,
				temperature;
				tol,
			)

			uk = @view eigen_k.vectors[:, n]
			ukq = @view eigen_kq.vectors[:, m]

			for alpha in 1:4, beta in 1:4
				chi[alpha, beta] += kernel * conj(uk[alpha]) * ukq[alpha] * conj(ukq[beta]) * uk[beta]
			end
		end
	end

	return chi / Nk^2
end

function diagonalized_kmesh(Nk::Integer, U)
	Nk > 0 || throw(ArgumentError("Nk must be positive"))

	energies = zeros(Float64, Nk, Nk, 4)
	vectors = zeros(ComplexF64, Nk, Nk, 4, 4)

	for ix in 1:Nk, iy in 1:Nk
		kx = 2pi * (ix - 1) / Nk
		ky = 2pi * (iy - 1) / Nk
		eigen_k = eigen(Hermitian(tight_binding(kx, ky, U, 0.0)))
		energies[ix, iy, :] = eigen_k.values
		vectors[ix, iy, :, :] = eigen_k.vectors
	end

	return (; Nk, energies, vectors)
end

function chi_matrix_q_index(qx_index::Integer, qy_index::Integer, cache, mu; temperature = 0.001, tol = 1e-8)
	Nk = cache.Nk
	chi = zeros(ComplexF64, 4, 4)

	for ix in 1:Nk, iy in 1:Nk
		ixq = mod1(ix + qx_index, Nk)
		iyq = mod1(iy + qy_index, Nk)

		for n in 1:4, m in 1:4
			kernel = susceptibility_kernel(
				cache.energies[ix, iy, n],
				cache.energies[ixq, iyq, m],
				mu,
				temperature;
				tol,
			)

			uk = @view cache.vectors[ix, iy, :, n]
			ukq = @view cache.vectors[ixq, iyq, :, m]

			for alpha in 1:4, beta in 1:4
				chi[alpha, beta] += kernel * conj(uk[alpha]) * ukq[alpha] * conj(ukq[beta]) * uk[beta]
			end
		end
	end

	return chi / Nk^2
end

function chi_matrix_R(R::Tuple{<:Integer, <:Integer}; Nk::Integer = 40, Nq::Integer = Nk, U = 0.0, mu = 0.0, temperature = 0.001, tol = 1e-8)
	Nq == Nk || throw(ArgumentError("Nq must equal Nk in this cached grid implementation"))

	cache = diagonalized_kmesh(Nk, U)
	chi_R = zeros(ComplexF64, 4, 4)

	for iqx in 0:(Nq-1), iqy in 0:(Nq-1)
		qx = 2pi * iqx / Nq
		qy = 2pi * iqy / Nq
		phase = exp(im * (qx * R[1] + qy * R[2]))
		chi_R += phase * chi_matrix_q_index(iqx, iqy, cache, mu; temperature, tol)
	end

	return chi_R / Nq^2
end
