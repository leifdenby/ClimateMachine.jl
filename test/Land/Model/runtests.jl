using Test, Pkg
using ClimateMachine.Land

@testset "Land" begin
    include("constant_hydraulic_conductivity.jl")
    include("constant_energy_model.jl")
end
