using JLD2
using Test

plot_script = joinpath(
	@__DIR__,
	"..",
	"scripts",
	"plotting",
	"plot_total_rkky.jl",
)

@testset "plot total RKKY" begin
	@test isfile(plot_script)

	if isfile(plot_script)
		include(plot_script)

		tmpdir = mktempdir()
		filepath = joinpath(tmpdir, "total.jld2")
		U_values = [0.0, 0.1, 0.2]
		rslt = (;
			U_values,
			J0 = [1.0, 1.1, 1.2],
			J1 = [0.0, 0.1, 0.2],
			J2 = [-1.0, -0.9, -0.8],
			J3 = [0.2, 0.1, 0.0],
		)
		jldsave(filepath; rslt)

		plots = plot_total_rkky_from_jld2(filepath; output_dir = tmpdir)

		@test sort(collect(keys(plots))) == ["J0", "J1", "J2", "J3"]
		for name in ("J0", "J1", "J2", "J3")
			outfile = joinpath(tmpdir, "rkky_total_$(name)_U0-8.png")
			@test isfile(outfile)
			@test filesize(outfile) > 0
		end
	end
end
