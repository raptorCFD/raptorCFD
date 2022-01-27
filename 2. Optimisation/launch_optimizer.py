from STEEPEST import opt_method_steepest
from PSO import opt_method_pso
from preliminary_study import preliminary_study
from input_check_opt import input_check_opt
from main_opt import *

input_check_opt(preliminary, preliminary_method, flow, n, objective_function, optim_method, typ,
                design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, TE,
                start_dist, start_alpha, start_crel, start_params, bounds_alpha, bounds_dist, bounds_crel,
                bounds_params, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo,
                external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining,
                ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size)

varbound = []
k = 0
for i in range(0, len(design_variables_alpha)):
    if design_variables_alpha[i] == 'Y':
        varbound = varbound + [bounds_alpha[k]]
        k = k + 1

for i in range(0, len(design_variables_dist)):
    k = 0
    for j in range(0, len(design_variables_dist[i])):
        if design_variables_dist[i][j] == 'Y':
            varbound = varbound + [bounds_dist[i][k]]
            k = k + 1
k = 0
for i in range(0, len(design_variables_crel)):
    if design_variables_crel[i] == 'Y':
        varbound = varbound + [bounds_crel[k]]
        k = k + 1

for i in range(0, len(design_variables_params)):
    k = 0
    for j in range(0, len(design_variables_params[i])):
        if design_variables_params[i][j] == 'Y':
            bounds_params[i][k]
            varbound = varbound + [bounds_params[i][k]]
            k = k + 1

if preliminary == 'YES':
    start_alpha, start_dist, start_crel, start_params, varbound = preliminary_study(preliminary_method, flow, n, objective_function, optim_method, typ, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel, start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size, varbound)

# Optimization
if optim_method == 'STEEPEST':
    opt_method_steepest(flow, n, objective_function, typ, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel, start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size, varbound, preliminary_method)

elif optim_method == 'PSO':
    opt_method_pso(flow, n, objective_function, typ, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel, start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size, varbound)
