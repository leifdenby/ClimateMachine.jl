using Test, Pkg
using ClimateMachine.Land

@testset "Land" begin
    include("test_hydraulic_conductivity_factors.jl")
    include("constant_moisture_model.jl")
    include("constant_energy_model.jl")
end
