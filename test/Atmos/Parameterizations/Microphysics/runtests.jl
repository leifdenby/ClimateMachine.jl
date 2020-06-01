using Test

include(joinpath(@__DIR__, "..","..","..", "testhelpers.jl"))

@testset "Microphysics tests" begin
    tests = [
        (1, "1_unit_tests.jl")
        (1, "2_saturation_adjustment.jl")
        (1, "3_warm_rain.jl")
    ]

    runmpi(tests, @__FILE__)
end
