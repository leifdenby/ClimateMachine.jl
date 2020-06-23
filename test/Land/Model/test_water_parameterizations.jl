using MPI
using OrderedCollections
using StaticArrays
using Test

using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using ClimateMachine
using ClimateMachine.Land.SoilWaterParameterizations

FT = Float64;

@test viscosity_factor(ConstantViscosity{FT}()) == 1
viscosity_model = TemperatureDependentViscosity{FT}(; T_ref = FT(1.0))
@test viscosity_factor(viscosity_model, FT(1.0)) == 1

@test impedance_factor(NoImpedance{FT}()) == 1
impedance_model = IceImpedance{FT}(; Ω = 2.0)
@test impedance_factor(impedance_model, 0.2, 0.4) == FT(0.1)

vg_model = vanGenuchten{FT}()
mm = MoistureDependent{FT}()
@test moisture_factor(mm, vg_model, FT(1)) == 1
@test moisture_factor(mm, vg_model, FT(0)) == 0

bc_model = BrooksCorey{FT}()
@test moisture_factor(mm, bc_model, FT(1)) == 1
@test moisture_factor(mm, bc_model, FT(0)) == 0

hk_model = Haverkamp{FT}()
@test moisture_factor(mm, hk_model, FT(1), FT(1)) == 1
@test moisture_factor(mm, hk_model, FT(0), FT(0)) == 1

@test moisture_factor(MoistureIndependent{FT}()) == 1

#Test hydraulic models of interest
#Constant K
@test hydraulic_conductivity(
    NoImpedance{FT}(),
    ConstantViscosity{FT}(),
    MoistureIndependent{FT}(),
) == 1

#Liquid model - haverkamp
@test hydraulic_conductivity(
    NoImpedance{FT}(),
    ConstantViscosity{FT}(),
    MoistureDependent{FT}(),
    Haverkamp{FT}(; A = FT(1.0));
    S_l = FT(0.5),
    ψ = FT(1.0),
) == FT(0.5)

#Liquid model - vG
@test hydraulic_conductivity(
    NoImpedance{FT}(),
    ConstantViscosity{FT}(),
    MoistureDependent{FT}(),
    vanGenuchten{FT}();
    S_l = FT(1.0),
) == 1

#Liquid and ice  - vG
@test hydraulic_conductivity(
    impedance_model,
    ConstantViscosity{FT}(),
    MoistureDependent{FT}(),
    vanGenuchten{FT}();
    θ_ice = 0.2,
    porosity = 0.4,
    S_l = 1.0,
) == FT(0.1)
#Liquid+Viscosity model - vG
@test hydraulic_conductivity(
    NoImpedance{FT}(),
    viscosity_model,
    MoistureDependent{FT}(),
    vanGenuchten{FT}();
    T = 1.0,
    S_l = 1.0,
) == 1
#Full model - vG
@test hydraulic_conductivity(
    impedance_model,
    viscosity_model,
    MoistureDependent{FT}(),
    vanGenuchten{FT}();
    θ_ice = 0.5,
    porosity = 1.0,
    T = 1.0,
    S_l = 1.0,
) == FT(0.1)

@test effective_saturation(0.5, -1.0) == 0
@test effective_saturation(0.5, 0.25) == 0.5

test_array = [0.5, 1.0]
n = FT(1.43)
m = 1.0 - 1.0 / n
α = FT(2.6)

@test pressure_head.(Ref(vg_model), Ref(1.0), Ref(0.001), test_array) ≈
      .-((-1 .+ test_array .^ (-1 / m)) .* α^(-n)) .^ (1 / n)
#test branching in pressure head
@test pressure_head(vg_model, 1.0, 0.001, 1.5) == 500

@test pressure_head.(Ref(hk_model), Ref(1.0), Ref(0.001), test_array) ≈
      .-((-1 .+ test_array .^ (-1 / m)) .* α^(-n)) .^ (1 / n)


m = FT(0.5)
ψb = FT(0.1656)
@test pressure_head(bc_model, 1.0, 0.001, 0.5) ≈ -ψb * 0.5^(-1 / m)
