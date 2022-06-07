"""
File with the constants for the Duns & Ros correlation
"""
from scipy.interpolate import interp1d


# Factors that describes limits of 1st Region (Bubble Flow + Slug Flow)
L1_x = [0, 16, 20, 30, 40, 50, 70, 275]
L1_y = [2, 2, 1.9, 1.6, 1.25, 1.1, 1, 1]

L2_x = [7.5, 20, 30, 40, 50, 60, 100, 275]
L2_y = [0.5, 0.75, 0.9, 1, 1.1, 1.1, 1.11, 1.11]


# Factors that needed to define liquid holdup
N_l_x = [0.002, 0.006, 0.007, 0.01, 0.02, 0.04, 0.05, 0.1, 0.2, 0.4, 2]
F1_y = [1.3, 1.3, 1.3, 1.3, 1.3, 1.5, 1.65, 2, 2, 1.8, 0.9]
F2_y = [0.25, 0.25, 0.25, 0.25, 0.28, 0.45, 0.6, 0.95, 1, 0.98, 0.7]
F3_y = [0.8, 0.9, 1, 1.3, 1.9, 2.5, 3, 3.25, 3.5, 3.75, 4]
F4_y = [-20, 5, 10, 25, 38, 50, 52, 55, 55, 55, 55]
F5_y = [0.22, 0.2, 0.19, 0.18, 0.17, 0.15, 0.13, 0.065, 0.048, 0.06, 0.11]
F6_y = [0.8, 0.05, 0, -0.1, -0.1, 0.65, 1, 2.1, 1.9, 1.8, 1.75]
F7_y = [0.14, 0.101, 0.099, 0.09, 0.07, 0.056, 0.052, 0.04, 0.034, 0.03, 0.025]

# Factors that needed to define friction factor for the liquid
coef_x = [0.001, 0.4, 0.7, 1, 2, 3, 6, 10, 20, 40, 100]
f2_y = [1, 1, 0.9, 0.75, 0.6, 0.5, 0.4, 0.35, 0.3, 0.25, 0.2]

# Interpolation by scipy
f1_func = interp1d(
    x=N_l_x,
    y=F1_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f2_func = interp1d(
    x=N_l_x,
    y=F2_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f3_func = interp1d(
    x=N_l_x,
    y=F3_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f4_func = interp1d(
    x=N_l_x,
    y=F4_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f5_func = interp1d(
    x=N_l_x,
    y=F5_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f6_func = interp1d(
    x=N_l_x,
    y=F6_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f7_func = interp1d(
    x=N_l_x,
    y=F7_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f2_y_func = interp1d(
    x=coef_x,
    y=f2_y,
    fill_value="extrapolate",
    kind="linear",
)
l1_func = interp1d(
    x=L1_x,
    y=L1_y,
    fill_value="extrapolate",
    kind="quadratic",
)
l2_func = interp1d(
    x=L2_x,
    y=L2_y,
    fill_value="extrapolate",
    kind="quadratic",
)

dr_const = {
    "f1": f1_func,
    "f2": f2_func,
    "f3": f3_func,
    "f4": f4_func,
    "f5": f5_func,
    "f6": f6_func,
    "f7": f7_func,
    "f2_y": f2_y_func,
    "l1": l1_func,
    "l2": l2_func,
}