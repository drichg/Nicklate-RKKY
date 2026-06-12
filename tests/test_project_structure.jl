using Test

project_root = normpath(joinpath(@__DIR__, ".."))
module_path = joinpath(project_root, "src", "NickelateRKKY.jl")

@testset "project structure" begin
	@test isfile(module_path)
	@test isdir(joinpath(project_root, "src", "models"))
	@test isdir(joinpath(project_root, "src", "susceptibility"))
	@test isdir(joinpath(project_root, "src", "ground_state"))
	@test isdir(joinpath(project_root, "src", "excitations"))

	if isfile(module_path)
		include(module_path)
		using .NickelateRKKY

		@test NickelateRKKY.PROJECT_ROOT == project_root
		@test NickelateRKKY.MODEL_DATA_DIR == joinpath(project_root, "data", "model")
		@test isfile(NickelateRKKY.Z_DATA_PATH)
		@test isfinite(NickelateRKKY.chemical_potential(0.0))
	end
end
