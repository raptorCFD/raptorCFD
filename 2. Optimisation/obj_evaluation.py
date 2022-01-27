import pandas as pd




def obj_evaluation_GA(X, matrix, n, objective_function, optim_method, typ, params, design_variables_alpha, design_variables_dist, crel, TE, start_dist, start_alpha, varbound, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, limit_exclusions):
    print('variables', X)
    params = [None] * nairfoils
    k = 0
    for i in range(0, len(params)):
        params[i] = [None] * len(bound_params[i])
        for j in range(0, len(bound_params[i])):
            if design_variables[i][j] == 'Y' and k != (len(X)):
                params[i][j] = X[k]
                k = k + 1
            else:
                params[i][j] = start_params[i][j]

    print('IGP Parameters', params)

    # wing_launch_shape_opt(n, params, alpha, dist, crel, typ, TE, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil)

    ITER = int(main_configuration(nairfoils, Re, Mach, mu, temp, alpha))

    # subp = subprocess.call("./CFD_run.sh")

    if objective_function[1] == 'LIFT':
        column_name = '      "CL"      '
    elif objective_function[1] == 'DRAG':
        column_name = '      "CD"      '
    else:
        column_name = '      "CEff"      '
    df = pd.read_csv("history.csv")

    obj_old = df.at[(df.index.stop - 1), column_name]
    sim_iter = df.index.stop - 1
    if column_name == '      "CL"      ':
        print('Last configuration has Lift Coefficient equal to: ' + str(obj_old) + '.\n')
    elif column_name == '      "CD"      ':
        print('Last configuration has Drag Coefficient equal to: ' + str(obj_old) + '.\n')
    else:
        print('Last configuration has Efficiency equal to: ' + str(obj_old) + '.\n')

    if objective_function[0] == 'MAX':
        obj = -obj_old
    else:
        obj = obj_old
        pass
    matrix.append([X, obj, sim_iter])

    if sim_iter == ITER:
        obj = 100000000

    return obj

