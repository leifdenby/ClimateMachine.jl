using Test, Pkg
@testset "Land" begin
    include("test_water_parameterizations.jl")
    include("constant_moisture_model.jl")
    include("constant_energy_model.jl")
end
