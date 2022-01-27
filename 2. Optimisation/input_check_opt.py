import sys
import numpy as np
from functools import partial


def _is_feasible_wrapper(func, x):
    return np.all(func(x) >= 0)


def _cons_none_wrapper(x):
    return np.array([0])


def _cons_ieqcons_wrapper(ieqcons, args, kwargs, x):
    return np.array([y(x, *args, **kwargs) for y in ieqcons])


def _cons_f_ieqcons_wrapper(f_ieqcons, args, kwargs, x):
    return np.array(f_ieqcons(x, *args, **kwargs))


def findBestFitProfile(x_target, y_target, TE, k, n, typ):
    import os
    import matplotlib.pyplot as plt
    from airfoilShape import airfoilIGP
    ieqcons = []
    f_ieqcons = None
    args = ()
    kwargs = {}
    S = 200 # 170     # swarm size
    D = 8  # the number of dimensions each particle has
    IGPpoints = 300 # 250  # Number of geometrical points for IGP shape tests' definition
    omega = 0.25
    phip = 0.5
    phig = 0.75
    maxiter = 20 #25
    minstep = 1e-4
    minfunc = 1e-4
    print('Searching (max %d iterations)...\n' % maxiter)
    lb = np.array([0.01, 0.02, -0.07, -0.10, 0.1, 0.02, 0.16, 0.14])
    ub = np.array([0.96, 0.97, 0.25, 0.21, 0.5, 0.33, 1.50, 4.90])
    vhigh = np.abs(ub - lb)
    vlow = -vhigh

    if f_ieqcons is None:
        if not len(ieqcons):
            cons = _cons_none_wrapper
        else:
            cons = partial(_cons_ieqcons_wrapper, ieqcons, args, kwargs)
    else:
        cons = partial(_cons_f_ieqcons_wrapper, f_ieqcons, args, kwargs)
    is_feasible = partial(_is_feasible_wrapper, cons)

    var = np.random.rand(S, D) * (ub - lb) + lb  # Initialize the particle's position

    v = np.zeros_like(var)  # particle velocities
    p = np.zeros_like(var)  # best particle positions
    fx = np.zeros(S)  # current particle function values
    fs = np.zeros(S, dtype=bool)  # feasibility of each particle
    fp = np.ones(S) * np.inf  # best particle function values
    g = []  # best swarm position
    fg = np.inf  # best swarm position starting value

    if not os.path.exists('Element_%d' % k):
        os.makedirs('Element_%d' % k)
    if typ == ('MY_FILE_%d' % (k + 1)):
        import shutil
        shutil.copy(('MY_FILE_%d.dat' % (k + 1)), 'Element_%d' % k)
    os.chdir('Element_%d' % k)

    data_id = open('IGP_shape_history.txt', 'w+')
    for i in range(0, len(var)):
        fx[i] = getErr(var[i], x_target, y_target, TE, IGPpoints)
        data_id.write('Candidate %d: [%.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f]. Fitness: %.8f.\n' % (
            i, var[i][0], var[i][1], var[i][2], var[i][3], var[i][4], var[i][5], var[i][6], var[i][7], fx[i]))
        fs[i] = is_feasible(var[i, :])

    # Store particle's best position (if constraints are satisfied)
    i_update = np.logical_and((fx < fp), fs)
    p[i_update, :] = var[i_update, :].copy()
    fp[i_update] = fx[i_update]

    # Update swarm's best position
    i_min = np.argmin(fp)
    if fp[i_min] < fg:
        fg = fp[i_min]
        g = p[i_min, :].copy()
    else:
        # At the start, there may not be any feasible starting point, so just
        # give it a temporary "best" point since it's likely to change
        g = var[0, :].copy()
    print('Iteration 0: best IGP shape is [%.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f].\n' % (g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]))
    data_id.write('Iteration 0: best IGP shape is [%.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f].\n' % (g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]))
    data_id.write('Iteration 0: objective function %.8f.\n' % fg)
    data_id.close()
    # Initialize the particle's velocity
    v = vlow + np.random.rand(S, D) * (vhigh - vlow)

    # Iterate until termination criterion met ##################################
    it = 1

    while it <= maxiter:
        rp = np.random.uniform(size=(S, D))
        rg = np.random.uniform(size=(S, D))

        # Update the particles velocities
        v = omega * v + phip * rp * (p - var) + phig * rg * (g - var)

        # Update the particles' positions
        var = var + v

        # Correct for bound violations
        maskl = var < lb
        masku = var > ub
        var = var * (~np.logical_or(maskl, masku)) + lb * maskl + ub * masku

        data_id = open('IGP_shape_history.txt', 'a')
        # Update objectives and constraints
        for i in range(0, len(var)):
            fx[i] = getErr(var[i], x_target, y_target, TE, IGPpoints)
            data_id.write('Candidate %d: [%.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f]. Fitness: %.8f.\n' % (
            i, var[i][0], var[i][1], var[i][2], var[i][3], var[i][4], var[i][5], var[i][6], var[i][7], fx[i]))
            fs[i] = is_feasible(var[i, :])

        # Store particle's best position (if constraints are satisfied)
        i_update = np.logical_and((fx < fp), fs)
        p[i_update, :] = var[i_update, :].copy()
        fp[i_update] = fx[i_update]

        # Compare swarm's best position with global best position
        i_min = np.argmin(fp)

        if fp[i_min] < fg:
            p_min = var[i_min, :].copy()
            stepsize = np.sqrt(np.sum((g - p_min) ** 2))

            if np.abs(fg - fp[i_min]) <= minfunc or stepsize <= minstep:
                print('Iteration %d: best IGP shape is [%.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f].\n' % (it,
                g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]))
                data_id.write('Iteration %d: best IGP shape is [%.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f].\n' % (it,
                g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]))
                data_id.write('Iteration %d: objective function %.8f.\n' % (it, fg))
            else:
                g = p_min.copy()
                fg = fx[i_min]
                print('Iteration %d: best IGP shape is [%.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f].\n' % (it,
                g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]))
                data_id.write('Iteration %d: best IGP shape is [%.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f].\n' % (it,
                                                                                                               g[0],
                                                                                                               g[1],
                                                                                                               g[2],
                                                                                                               g[3],
                                                                                                               g[4],
                                                                                                               g[5],
                                                                                                               g[6],
                                                                                                               g[7]))

                data_id.write('Iteration %d: objective function %.8f.\n' % (it, fg))
        else:
            data_id.write('Iteration %d: best IGP shape is [%.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f].\n' % (it,
                                                                                                           g[0],
                                                                                                           g[1],
                                                                                                           g[2],
                                                                                                           g[3],
                                                                                                           g[4],
                                                                                                           g[5],
                                                                                                           g[6],
                                                                                                           g[7]))
            data_id.write('Iteration %d: objective function %.8f.\n' % (it, fg))
        it += 1
        data_id.close()
    xcurr, ycurr, TE[k] = airfoilIGP(g, n, TE[k])  # x,y parameterized by IGP

    plt.figure()
    plt.plot(x_target, y_target, 'k', linewidth=3, label='Input shape')
    plt.plot(xcurr, ycurr, 'b', linewidth=3, label='IGP shape')
    plt.axis('equal')
    plt.legend()
    plt.savefig('airfoil_%d_IGP_definition.png' % k)
    plt.close()
    path_parent = os.path.dirname(os.getcwd())
    os.chdir(path_parent)

    return g


def getErr(p, x_target, y_target, TE, IGPpoints):
    import math as mt
    import numpy as np
    from airfoilShape import airfoilIGP
    x, y, TE = airfoilIGP(p, IGPpoints, TE)

    J = 0

    #if n < IGPpoints:
    for ii in range(0, len(x_target)):
        idx = np.argmin([mt.sqrt((x_target[ii] - x[jj])**2 + (y_target[ii] - y[jj])**2) for jj in range(0, len(x))])
        J = J + mt.sqrt((x[idx] - x_target[ii])**2 + (y[idx] - y_target[ii])**2)

    #else:
        #for ii in range(0, len(x)):
        #    idx = np.argmin([mt.sqrt((x_target[jj] - x[ii]) ** 2 + (y_target[jj] - y[ii]) ** 2) for jj in range(0, len(x_target))])
        #    J = J + mt.sqrt((x[ii] - x_target[idx]) ** 2 + (y[ii] - y_target[idx]) ** 2)

    return J


def input_check_opt(preliminary, preliminary_method, flow, n, objective_function, optim_method, typ, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel, start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size):

    print('Checking consistency of given input values...\n')

    # Chosen parametrization
    for i in range(0, len(typ)):
        if typ[i] != 'IGP' and typ[i] != 'NACA4' and typ[i] != 'NACA5' and typ[i] != 'NACA4_MOD' and typ[
            i] != 'NACA16' and typ[i] != 'BENZING' and typ[i] != ('MY_FILE_%d' % (i + 1)):
            print(
                "Error: the chosen parametrization should be one of the following: IGP, NACA4, NACA5, NACA4_MOD, NACA16 or MY_FILE_%d for airfoil n° %d. \n The written should be as shown above. Check the item: 'typ'." % (
                i + 1, i + 1))
            sys.exit()
            # Option "BICONVEX" to be developed. "BENZING" to be completed

    if n <= 0:
        print(
            "Error: number of points for geometry generation should be a value greater than 0. Hints: 200 or 500. \nCheck the following item: 'n'.")
        sys.exit()

    for i in range(0, 2):
        if i == 0:
            if objective_function[i] != 'MAX' and objective_function[i] != 'MIN':
                print("Error: first term objective function item defines maximum or minimum objective function typology; then, fist term accepts only string 'MAX' or 'MIN'. \nCheck the following item: 'objective_function'.")
                sys.exit()
        else:
            if objective_function[i] != 'LIFT' and objective_function[i] != 'DRAG' and objective_function[i] != 'EFFICIENCY':
                print("Error: from second term to last of objective function item define objective function typology; ; then, fist term accepts only string 'LIFT' or 'DRAG' or 'EFFICIENCY'. \nCheck the following item: 'objective_function'.")
                sys.exit()




    # Length of each parameter
    if len(objective_function) != 2:
        print("Error: objective function item accepts only two terms as input, which are max/min and type of objective function; then, length of this item must be 2. \nCheck the following item length: 'objective_function'.")
        sys.exit()
    elif len(design_variables_alpha) != len(typ):
        print(
            "Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \nCheck the following item length: 'design_variables_alpha'.")
        sys.exit()
    elif len(design_variables_dist) != (len(typ) - 1):
        print("Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \nCheck the following item length: 'design_variables_dist'.")
        sys.exit()
    elif len(design_variables_crel) != (len(typ) - 1):
        print("Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \nCheck the following item length: 'design_variables_crel'.")
        sys.exit()
    elif len(design_variables_params) != len(typ):
        print("Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \nCheck the following item length: 'design_variables_params'.")
        sys.exit()
    elif len(start_alpha) != len(typ):
        print(
            "Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \n Check the following item length: 'start_alpha'.")
        sys.exit()
    elif len(start_dist) != len(typ):
        print(
            "Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \n Check the following item length: 'start_dist'.")
        sys.exit()
    elif len(start_crel) != len(typ):
        print(
            "Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \n Check the following item length: 'start_crel'.")
        sys.exit()
    elif len(start_params) != len(typ):
        print("Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \nCheck the following item length: 'start_params'.")
        sys.exit()
    elif len(y_plus) != len(typ):
        print("Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \n Check the following item length: 'y_plus'.")
        sys.exit()
    elif len(thick) != len(typ):
        print("Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \n Check the following item length: 'thick'.")
        sys.exit()
    elif len(progr) != len(typ):
        print("Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \n Check the following item length: 'progr'.")
        sys.exit()
    elif len(wall_refining) != len(typ):
        print("Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \n Check the following item length: 'wall_refining'.")
        sys.exit()
    else:
        pass

    if optim_method != 'GA' and optim_method != 'PSO' and optim_method != 'STEEPEST':
        print("Error: optimization method item accepts only one string term between the followings: 'GA', 'PSO', 'STEEPEST'. \n Check the following item: 'optim_method'.")
        sys.exit()

    for i in range(0, len(start_params)):
        for j in range(0, len(start_params[i])):
            if type(start_params[i][j]) != int and type(start_params[i][j]) != float and start_params[i][j] != ('MY_FILE_%d' % (i+1)):
                print("Error: each value inserted in 'start_params' must be an integer or float value. \n Check the item: 'start_params'.")
                sys.exit()
            else:
                pass

        if design_variables_alpha[i] != 'Y' and design_variables_alpha[i] != 'N':
            print("Error: each value inserted in 'design_variables_alpha' must be defined by string (options for each term: 'Y', 'N'). \n Check the item: 'design_variables_alpha'.")
            sys.exit()
        else:
            pass

        if i != 0 and design_variables_crel[i-1] != 'Y' and design_variables_crel[i-1] != 'N':
            print("Error: each value inserted in 'design_variables_crel' must be defined by string (options for each term: 'Y', 'N'). \n Check the item: 'design_variables_crel'.")
            sys.exit()
        else:
            pass
        count = 0
        if typ[i] == "IGP":
            for j in range(0, len(bounds_params[count])):
                for k in range(0, len(bounds_params[count][j])):
                    if type(bounds_params[count][j][k]) != int and type(bounds_params[i][j][k]) != float:
                        print("Error: each value inserted in 'bounds_params' must be an integer or float value. \n Check the item: 'bounds_params'.")
                        sys.exit()


    for i in range(0, len(design_variables_dist)):
        for j in range(0, len(design_variables_dist[i])):
            if design_variables_dist[i][j] != 'Y' and design_variables_dist[i][j] != 'N':
                print("Error: each value inserted in 'design_variables_dist' must be defined by string (options for each term: 'Y', 'N'). \n Check the item: 'design_variables_dist'.")
                sys.exit()
            else:
                pass

    # Distance between airfoils item
    for i in range(0, len(typ)):
        if (type(start_alpha[i]) != int and type(start_alpha[i]) != float) or ((type(start_alpha[i]) == int or type(start_alpha[i]) == float) and (start_alpha[i] > 90 or start_alpha[i] < -90)):
            print(
                "Error: angle of attack of each airfoil should be written as integer or float number. \n Check the item: 'alpha'.")
            sys.exit()

        if (type(start_crel[i]) != int and type(start_crel[i]) != float) or start_crel[i] <= 0:
            print("Error: relative chord of each airfoil should be defined by an integer or float number greater than 0. \n Check the item: 'crel'.")
            sys.exit()

    for i in range(0, len(start_dist)):
        if len(start_dist[i]) != 2:
            print("Error: 'start_dist' term number (" + str(i) + ") must have length equal to 2.\n Check the item: 'start_dist'.")
            sys.exit()
        for j in range(0, len(start_dist[i])):
            if type(start_dist[i][j]) != int and type(start_dist[i][j]) != float:
                print(
                    "Error: distance between airfoils should be an integer or float number. \n Check the item: 'start_dist'.")
                sys.exit()

    # for i in range(0, len(bounds_dist)):
    #     if len(bounds_dist[i]) != 2:
    #         print("Error: 'bounds_dist' term number (" + str(i) + ") must have length equal to 2.\n Check the item: 'bounds_dist'.")
    #         sys.exit()
    #     for j in range(0, len(bounds_dist[i])):
    #         if type(bounds_dist[i][j]) != int and type(bounds_dist[i][j]) != float:
    #             print("Error: each value inserted in 'bounds_dist' must be an integer or float value. \n Check the item: 'bounds_dist'.")
    #             sys.exit()



    # Reynolds number and dynamic viscosity
    if Re <= 0:
        print("Error: Reynolds number should be greater than 0. \n Check the item: 'Re'.")
        sys.exit()
    elif mu <= 0:
        print("Error: dynamic viscosity should be greater than 0 [Pa * s]. \n Check the item: 'mu'.")
        sys.exit()
    elif type(Re) != int and type(Re) != float:
        print("Error: Reynolds number should be an integer or float. \n Check the following item: 'Re'.")
        sys.exit()
    elif type(mu) != int and type(mu) != float:
        print("Error: dynamic viscosity should be an integer or float. \n Check the following item: 'mu'.")
        sys.exit()
    else:
        pass

    # Nodes
    if nodes <= 0:
        print("Error: the number of airfoils must be greater than 0.")
        sys.exit()
    elif type(nodes) != int and type(nodes) != float:
        print("Error: the number of airfoils must be an integer or a float number.")
        sys.exit()
    else:
        pass

    for i in range(0, len(y_plus)):
        if (type(y_plus[i]) != int and type(y_plus[i]) != float) or (y_plus[i] <= 0 or y_plus[i] > 1):
            print('Error: dimensionless wall distance "y_plus" must be a vector of numerical values, greater than 0 and lower than 1. \n Check item "y_plus".')
            sys.exit()
        else:
            pass


    for i in range(0, len(thick)):
        if (type(thick[i]) != int and type(thick[i]) != float) or thick[i] <= 0:
            print('Error: structured region thickness must be a vector of numerical values, each greater than 0. Pay attention to avoid intersections between structured regions and airfoils\' geometries. \n Check item "thick".')
            sys.exit()
        else:
            pass

    for i in range(0, len(progr)):
        if ((type(progr[i]) != int and type(progr[i]) != float)) or (progr[i] < 1 or progr[i] > 1.3):
            print('Error: structured region progression item must be a vector of numerical values, each one greater than 1 and lower than 1.3. \n Check item "progr".')
            sys.exit()
        else:
            pass


    if (type(temp) != int and type(temp) != float):
        print('Error: free-stream temperature item must be a numerical value. \n Check item "temp".')
        sys.exit()
    else:
        pass

    if (type(Mach) != int and type(Mach) != float and Mach <= 0):
        print('Error: free-stream Mach item must be a numerical value greater than 0. \n Check item "Mach".')
        sys.exit()
    else:
        pass


    if type(ellipse_dimension) != float and ellipse_dimension < 0.1 and ellipse_dimension > 0.95:
        print('Error: dimension of the ellipse should be defined by a float number between 0.1 and 0.95. \n Check item "ellipse_dimension.')
        sys.exit()


    if (type(Mesh_Algo) != int and Mesh_Algo != 'DEFAULT') or (type(Mesh_Algo) == int and (Mesh_Algo <= 0 or Mesh_Algo >= 10)):
        print('Error: mesh algorithm item accepts only integer values between 1 and 9, or the string "DEFAULT" which corresponds to integer 6 (Frontal-Delaunay). \n Check item "Mesh_Algo".')
        sys.exit()

    if external_pts != 'YES' and external_pts != 'NO':
        print('Error: ellipse\'s external points\' item "external_pts" accepts only strings "YES" or "NO". \n Check item "external_pts".')
        sys.exit()

    if type(wake_length) != int and type(wake_length) != float and wake_length < 1:
        print('Error: wake length item should be an numerical value greater than 1. It is also suggested a value lower than 250-500 approximately (uncertain since it strictly depends on dimension of geometry chosen). \n Check item "wake_length".')
        sys.exit()

    if type(wake_progr) != int and type(wake_progr) != float and wake_progr < 1:
        print('Error: wake element progression item should be an numerical value greater than or equal to 1. \n Check item "wake_progr".')
        sys.exit()

    if type(semicircle_dimension) != int and type(semicircle_dimension) != float and semicircle_dimension <= 1:
        print('Error: external semi-circle item should be an numerical value greater than 1. \n Check item "semicircle_dimension".')
        sys.exit()

    if type(semicircle_elem_factor) != int and type(semicircle_elem_factor) != float and semicircle_dimension <= 1:
        print('Error: external semi-circle item should be an numerical value greater than 1. \n Check item "semicircle_elem_factor".')
        sys.exit()

    for i in range(0, len(wall_refining)):
        if (type(wall_refining[i]) != int and type(wall_refining[i]) != float) or wall_refining[i] <= 0:
            print('Error: wall refining\'s item should be an numerical value greater than 0. \n Check item "wall_refining".')
            sys.exit()

    if (type(ellipse_refining) != int and type(ellipse_refining) != float) or ellipse_refining <= 0:
        print('Error: ellipse refining\'s item should be an numerical value greater than 0. \n Check item "ellipse_refining".')
        sys.exit()

    if type(ref_airfoil) != int and ref_airfoil != "DEFAULT" and ref_airfoil < 1 and ref_airfoil > len(typ):
        print(
            'Error: reference airfoil\'s item should be an integer value greater or equal to 1, or string written "DEFAULT". \n Check item "ref_airfoil".')
        sys.exit()

    if limit_exclusions[0] != 'YES' and limit_exclusions[0] != 'NO':
        print(
            'Error: item related to individuals which could reach upper limit of total iterations in RANS (and possibly exclude them) must be a string written "YES" or "NO". \n Check item "limit_exclusions".')
        sys.exit()
    if limit_exclusions[0] == 'YES':
        if limit_exclusions[1] >= 1:
            print(
                'Error: item related to individuals which could reach upper limit of total iterations in RANS (and possibly exclude them) in 2nd index must be float value lower than 1. \n Check item "limit_exclusions".')
            sys.exit()
        if dev_std <= 0:
            print(
                'Error: iter related to standard deviation to be evaluate exclusions in optimisation  (read only with "limit_exclusions" with "YES" at 1st index) must be a numerical value greater than 0. \n Check item "dev_std".')
            sys.exit()
        if iter_range <= 0:
            print(
                'Error: iter related to range of iteration to be evaluate exclusions in optimisation (read only with "limit_exclusions" with "YES" at 1st index) must be a numerical value greater than 0. \n Check item "iter_range".')
            sys.exit()

    if type(n) != int or (type(n) == int and n <= 0):
        print('Error: number of generating points for each airfoil should be an integer value greater than 0. \n Check item "n".')
        sys.exit()

    if preliminary != 'YES' and preliminary != 'NO':
        print('Error: item "preliminary" must be defined by string options "YES" or "NO". \n Check item "preliminary".')
        sys.exit()
    elif preliminary == 'YES' and preliminary_method[0] != 'EULER' and preliminary_method[0] != 'VISCOUS' and preliminary_method[0] != 'HS':
        print('Error: item "preliminary_method" at position 0 must be defined by string options "EULER" or "VISCOUS" or "HS". \n Check item "preliminary_method" at first list position.')
        sys.exit()
    elif preliminary == 'YES' and preliminary_method[1] != 'BB' and preliminary_method[1] != 'LHS':
        print('Error: item "preliminary_method" at position 1 must be defined by string options "BB" or "LHS". \n Check item "preliminary_method" at second list position.')
        sys.exit()

    if type(farfield_size) != int or (type(farfield_size) == int and farfield_size <= sum(start_crel)):
        print('Error: item "farfield_size" must be an integer value higher than sum of relative chords. \n Check item "farfield_size".')
        sys.exit()

    if flow != 'EULER' and flow != 'VISCOUS':
        print('Error: item "flow" must be defined by string options "EULER" or "VISCOUS". \n Check item "flow".')
        sys.exit()

    # IMPOSING IGP PARAMETRISATION WHENEVER SHAPE OPTIMISATION IS INVOLVED
    for k in range(0, len(typ)):
        for var in design_variables_params[k]:
            if var == 'Y':
                if typ[k] != "IGP":
                    print('IGP shape research for airfoil (%d).\n' % k)
                    if typ[k] == "NACA4":
                        from airfoilShape import airfoilNACA4
                        x_target, y_target = airfoilNACA4(start_params[k], n, TE[k])  # x,y parameterized using NACA 4-digit equation
                    elif typ[k] == "NACA5":
                        from airfoilShape import airfoilNACA5
                        x_target, y_target = airfoilNACA5(start_params[k], n, TE[k])  # x,y parameterized using NACA 5-digit equation
                    elif typ[k] == "NACA4_MOD":
                        from airfoilShape import airfoilNACA4_MOD
                        x_target, y_target = airfoilNACA4_MOD(start_params[k], n, TE[k])  # x,y parameterized using NACA 4-digit modified equation
                    elif typ[k] == "NACA16":
                        from airfoilShape import airfoilNACA16
                        x_target, y_target = airfoilNACA16(start_params[k], n, TE[k])
                    elif typ[k] == "BENZING":
                        from airfoilShape import airfoilBENZING
                        x_target, y_target, TE[k] = airfoilBENZING(start_params[k], n, TE[k], k)
                    elif typ[k] == ("MY_FILE_%d" % (k + 1)):
                        from airfoilShape import airfoilMY_FILE
                        x_target, y_target, TE[k] = airfoilMY_FILE(k + 1, TE[k])
                    start_params[k] = findBestFitProfile(x_target, y_target, TE, k, n, typ[k])
                    print('Found! Check folder "Element %d" for further details.\n' % k)
                    typ[k] = "IGP"
                break

    print('Input check completed! \n')

    return []
