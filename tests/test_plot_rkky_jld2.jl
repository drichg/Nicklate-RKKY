using JLD2
using Test

plot_script = joinpath(
	@__DIR__,
	"..",
	"scripts",
	"plotting",
	"plot_orbital_rkky.jl",
)

@testset "plot orbital RKKY" begin
	@test isfile(plot_script)

	if isfile(plot_script)
		include(plot_script)

		tmpdir = mktempdir()
		filepath = joinpath(tmpdir, "orbital.jld2")
		U_values = [0.0, 0.1, 0.2]
		rslt = (;
			U_values,
			J0 = [1.0 2.0 3.0; 1.1 2.1 3.1; 1.2 2.2 3.2],
			J1 = [2.0 4.0 6.0; 3.0 6.0 9.0; 4.0 8.0 12.0],
			J2 = [0.0 1.0 2.0; 0.1 1.1 2.1; 0.2 1.2 2.2],
			J3 = [1.0 2.0 3.0; 1.0 2.0 3.0; 2.0 4.0 6.0],
		)
		jldsave(filepath; rslt)

		cd(tmpdir) do
			plots, loaded = plot_rkky_from_jld2(filepath; output_dir = tmpdir)
			@test loaded.U_values == U_values
			@test sort(collect(keys(plots))) == ["J0", "J1", "J2", "J3", "J3/J1"]
			for name in ("J0", "J1", "J2", "J3")
				outfile = "rkky_$(name)_U0-8.png"
				@test isfile(outfile)
				@test filesize(outfile) > 0
			end
			@test isfile("rkky_J3_over_J1_U0-8.png")
			@test filesize("rkky_J3_over_J1_U0-8.png") > 0
		end
	end
end
