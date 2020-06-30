using Test, Pkg
@testset "Land" begin
    include("test_water_parameterizations.jl")
    include("constant_moisture_model.jl")
    include("haverkamp_test.jl")
end
