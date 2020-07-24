using ClimateMachine
const clima_dir = dirname(dirname(pathof(ClimateMachine)));
include(joinpath(clima_dir, "tutorials", "Numerics", "DGMethods", "Box1D.jl"))

run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_no_filter.pdf")

run_box1D(4, 0.0, 1.0, 1.0, "box_1D_4_tmar.pdf", tmar_filter = true)

run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_cutoff_1.pdf",
    cutoff_filter = true,
    cutoff_param = 1,
)

run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_cutoff_3.pdf",
    cutoff_filter = true,
    cutoff_param = 3,
)

run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_exp_1_4.pdf",
    exp_filter = true,
    exp_param_1 = 1,
    exp_param_2 = 4,
)

run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_exp_1_8.pdf",
    exp_filter = true,
    exp_param_1 = 1,
    exp_param_2 = 8,
)

run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_exp_1_32.pdf",
    exp_filter = true,
    exp_param_1 = 1,
    exp_param_2 = 32,
)

run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_boyd_1_4.pdf",
    boyd_filter = true,
    boyd_param_1 = 1,
    boyd_param_2 = 4,
)

run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_boyd_1_8.pdf",
    boyd_filter = true,
    boyd_param_1 = 1,
    boyd_param_2 = 8,
)

run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_boyd_1_32.pdf",
    boyd_filter = true,
    boyd_param_1 = 1,
    boyd_param_2 = 32,
)

run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_tmar_exp_1_8.pdf",
    exp_filter = true,
    tmar_filter = true,
    exp_param_1 = 1,
    exp_param_2 = 8,
)

run_box1D(
    4,
    0.0,
    1.0,
    1.0,
    "box_1D_4_tmar_boyd_1_8.pdf",
    boyd_filter = true,
    tmar_filter = true,
    boyd_param_1 = 1,
    boyd_param_2 = 8,
)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

