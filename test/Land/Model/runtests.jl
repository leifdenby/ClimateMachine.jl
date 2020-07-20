using Test, Pkg
@testset "Land" begin
    include("test_water_parameterizations.jl")
    include("constant_moisture_model.jl")
<<<<<<< HEAD
    include("haverkamp_test.jl")
    include("test_bc.jl")#this maybe isnt necessary?
=======
>>>>>>> master
end
