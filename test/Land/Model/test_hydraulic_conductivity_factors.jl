using MPI
using OrderedCollections
using StaticArrays
using CLIMAParameters
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
using ClimateMachine
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.DGMethods: BalanceLaw, LocalGeometry
using ClimateMachine.MPIStateArrays
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.VariableTemplates
using ClimateMachine.SingleStackUtils
using Test
FT = Float64;
ClimateMachine.init(; disable_gpu = true);

const clima_dir = dirname(dirname(pathof(ClimateMachine)));
include(joinpath(clima_dir, "src", "Land","Model","soil_water.jl"));

@test viscosity_factor(ConstantViscosity{FT}()) == 1
viscosity_model = TemperatureDependentViscosity{FT}(; T_ref = FT(1.0))
@test viscosity_factor(viscosity_model, FT(1.0)) == 1

@test impedance_factor(NoImpedance{FT}()) == 1
impedance_model = IceImpedance{FT}(;Ω = 2.0)
@test impedance_factor(impedance_model, 0.2, 0.4) == FT(0.1)

vg_model = vanGenuchten{FT}()
mm = MoistureDependent{FT}()
@test moisture_factor(mm, vg_model, FT(1)) == 1
@test moisture_factor(mm, vg_model, FT(0)) == 0

bc_model = BrooksCorey{FT}()
@test moisture_factor(mm, bc_model,FT(1)) == 1
@test moisture_factor(mm, bc_model,FT(0)) == 0

hk_model = Haverkamp{FT}()
@test moisture_factor(mm, hk_model,FT(1),FT(1)) == 1
@test moisture_factor(mm, hk_model,FT(0),FT(0)) == 1

@test moisture_factor(MoistureIndependent{FT}()) == 1

#Test hydraulic models of interest
#Constant K
@test hydraulic_conductivity(NoImpedance{FT}(), ConstantViscosity{FT}(), MoistureIndependent{FT}();) == 1

#Liquid model - haverkamp
@test hydraulic_conductivity(NoImpedance{FT}(), ConstantViscosity{FT}(), MoistureDependent{FT}(), Haverkamp{FT}(; A = FT(1.0)); S_l =  FT(0.5), ψ = FT(1.0)) == FT(0.5)

#Liquid model - vG
@test hydraulic_conductivity(NoImpedance{FT}(), ConstantViscosity{FT}(), MoistureDependent{FT}(), vanGenuchten{FT}(); S_l = FT(1.0)) == 1

#Liquid and ice  - vG
@test hydraulic_conductivity(impedance_model, ConstantViscosity{FT}(), MoistureDependent{FT}(),  vanGenuchten{FT}();  θ_ice = 0.2, porosity = 0.4, S_l = 1.0) == FT(0.1)
#Liquid+Viscosity model - vG
@test hydraulic_conductivity(NoImpedance{FT}(), viscosity_model, MoistureDependent{FT}(),  vanGenuchten{FT}(); T = 1.0, S_l = 1.0) == 1
#Full model - vG
@test hydraulic_conductivity(impedance_model, viscosity_model, MoistureDependent{FT}(), vanGenuchten{FT}();  θ_ice = 0.5, porosity = 1.0, T = 1.0, S_l = 1.0 ) == FT(0.1)

