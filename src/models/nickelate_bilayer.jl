function z_data_row(U)
	return Int(round(U / 0.1 + 1))
end

function tight_binding(kx, ky, U, mu)
	t1 = 0
	t2 = 0
	t1z = 0
	t2z = 0
	t11x = -0.1123
	t22x = -0.4897
	t12x = 0.2425
	t11xy = -0.0142
	t22xy = 0.0686
	t11xz = 0.0257
	t22xz = 0.0006
	t12xz = -0.0370

	h11 = t11x * 2 * (cos(kx) + cos(ky)) + t11xy * 2 * (cos(kx + ky) + cos(kx - ky))
	h22 = t22x * 2 * (cos(kx) + cos(ky)) + t22xy * 2 * (cos(kx + ky) + cos(kx - ky))
	h12 = t12x * 2 * (cos(kx) - cos(ky))
	h13 = t11xz * 2 * (cos(kx) + cos(ky))
	h24 = t22xz * 2 * (cos(kx) + cos(ky))
	h14 = t12xz * 2 * (cos(kx) - cos(ky))

	h1 = [
		h22 h12 h24 h14;
		h12 h11 h14 h13;
		h24 h14 h22 h12;
		h14 h13 h12 h11
	]
	h2 = [
		t2 0 t2z 0;
		0 t1 0 t1z;
		t2z 0 t2 0;
		0 t1z 0 t1
	]

	h = renormalize(h1, h2, U)
	h -= mu * I(4)
	return h
end

function renormalize(h1, h2, U)
	data = readdlm(Z_DATA_PATH, '\t')
	row = z_data_row(U)
	z = data[row, 2:5]
	onsite_shift = data[row, 6:9]

	Z = diagm(z)
	onsite_shift = diagm(onsite_shift)
	u = 1 / sqrt(2)
	basis_transform = [
		u u 0 0;
		0 0 u -u;
		u -u 0 0;
		0 0 u u
	]

	H1 = inv(basis_transform) * h1 * basis_transform
	H2 = inv(basis_transform) * h2 * basis_transform
	H1 = Z * H1 * Z
	H2 = onsite_shift
	H = basis_transform * H1 * inv(basis_transform) + basis_transform * H2 * inv(basis_transform)

	if !isapprox(H, adjoint(H))
		println("Error!!! non hermitian matrix is generated")
	end
	return H
end

function chemical_potential(U)
	data = readdlm(Z_DATA_PATH, '\t')
	return data[z_data_row(U), 10]
end
