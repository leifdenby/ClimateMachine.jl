using Test, Pkg
@testset "Land" begin
    #include("test_heat_parameterizations.jl")
    include("Bonan_temperature_test.jl")
    #include("haverkamp_test.jl")
    #include("test_bc.jl")
    #include("prescribed_twice.jl")
end
