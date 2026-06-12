using LinearAlgebra
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))
using DelimitedFiles
using Plots
include(joinpath(@__DIR__, "..", "..", "src", "NickelateRKKY.jl"))
using .NickelateRKKY
Nq=50
U = 0.0
cache = diagonalized_kmesh(Nq, U)
mu = find_chemical_potential(1.5, 4.0; Nk = Nq)
chi11_q = zeros(Float64, Nq, Nq)
chi22_q = zeros(Float64, Nq, Nq)

for iqx in 0:(Nq-1), iqy in 0:(Nq-1)
	chi = chi_matrix_q_index(iqx, iqy, cache, mu;)
	chi11_q[iqx+1, iqy+1] = real(chi[1, 1])
	chi22_q[iqx+1, iqy+1] = real(chi[2, 2])
end
qvals = 2π .* (0:(Nq-1)) ./ Nq
function fftshift2(A)
	N1, N2 = size(A)
	return circshift(A, (N1 ÷ 2, N2 ÷ 2))
end

qshift = 2π .* ((-Nq÷2):(Nq÷2-1)) ./ Nq

chi11_shift = fftshift2(chi11_q)
chi22_shift = fftshift2(chi22_q)

p1 = heatmap(
	qshift,
	qshift,
	permutedims(chi11_shift);
	xlabel = "qₓ",
	ylabel = "qᵧ",
	title = "Re χ₁₁(q), shifted BZ",
	aspect_ratio = 1,
	colorbar_title = "χ₁₁",
)

p2 = heatmap(
	qshift,
	qshift,
	permutedims(chi22_shift)+permutedims(chi11_shift);
	xlabel = "qₓ",
	ylabel = "qᵧ",
	title = "Re χ₂₂(q), shifted BZ",
	aspect_ratio = 1,
	colorbar_title = "χ₂₂",
)

plot(p1, p2; layout = (1, 2), size = (1000, 420))
# savefig("chi11_chi22_q_shifted_U$(U)_Nq$(Nq).png")
