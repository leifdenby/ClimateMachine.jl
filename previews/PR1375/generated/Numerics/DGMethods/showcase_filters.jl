using ClimateMachine
const clima_dir = dirname(dirname(pathof(ClimateMachine)));
include(joinpath(clima_dir, "tutorials", "Numerics", "DGMethods", "Box1D.jl"))

run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_no_filter.svg");

run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_tmar.svg", tmar_filter = true);

run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_cutoff_1.svg",
    cutoff_filter = true,
    cutoff_param = 1,
);

# `CutoffFilter(grid, Nc=3)`:
#run_box1D(

#)
# ![](box_1D_4_cutoff_3.svg)

# `ExponentialFilter(grid, Nc=1, s=4)`:
#run_box1D(

#)
# ![](box_1D_4_exp_1_4.svg)

# `ExponentialFilter(grid, Nc=1, s=8)`:
#run_box1D(

#)
# ![](box_1D_4_exp_1_8.svg)

# `ExponentialFilter(grid, Nc=1, s=32)`:
#run_box1D(

#)
# ![](box_1D_4_exp_1_32.svg)

# `BoydVandevenFilter(grid, Nc=1, s=4)`:
#run_box1D(

#)
# ![](box_1D_4_boyd_1_4.svg)

# `BoydVandevenFilter(grid, Nc=1, s=8)`:
#run_box1D(

#)
# ![](box_1D_4_boyd_1_8.svg)

# `BoydVandevenFilter(grid, Nc=1, s=32)`:
#run_box1D(

#)
# ![](box_1D_4_boyd_1_32.svg)

# `ExponentialFilter(grid, Nc=1, s=8)` and `TMARFilter()`:
#run_box1D(

#)
# ![](box_1D_4_tmar_exp_1_8.svg)

# `BoydVandevenFilter(grid, Nc=1, s=8)` and `TMARFilter()`:
#run_box1D(

#)
# ![](box_1D_4_tmar_boyd_1_8.svg)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

