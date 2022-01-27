import time
from random import randrange

import numpy as np
from cfg_printing import cfg_printing
from main_opt_method import main_PSO_algo_params
from main_configuration_file import *
import subprocess
import shutil
import pandas as pd
from wing_launch_opt import wing_launch_opt
import os
# from obj_evaluation import obj_evaluation_PSO
from functools import partial
import math as mt
import multiprocessing as mp
from lift_driver_mode import *


def _is_feasible_wrapper(func, x):
    return np.all(func(x) >= 0)


def _cons_none_wrapper(x):
    return np.array([0])


def _cons_ieqcons_wrapper(ieqcons, args, kwargs, x):
    return np.array([y(x, *args, **kwargs) for y in ieqcons])


def _cons_f_ieqcons_wrapper(f_ieqcons, args, kwargs, x):
    return np.array(f_ieqcons(x, *args, **kwargs))

def obj_evaluation_PSO(objective_function, limit_exclusions, iter_range, dev_std, j, x):
    os.chdir('./Process_%d' % j)
    path2 = ("history.csv")

    if objective_function[1] == 'LIFT':
        column_name = '       "CL"       '
    elif objective_function[1] == 'DRAG':
        column_name = '       "CD"       '
    else:
        column_name = '      "CEff"      '

    try:
        df = pd.read_csv(path2)
        obj_old = df.at[(df.index.stop - 1), column_name]
        sim_iter = df.index.stop - 1

        if sim_iter <= iter_range:
            iter_range = sim_iter - 1

        sum_mean = 0
        for i in range(sim_iter - iter_range, sim_iter):
            sum_mean = sum_mean + df.at[i, column_name]
        mean_obj = sum_mean / iter_range

        sum_dev_std = 0
        for i in range(sim_iter - iter_range, sim_iter):
            sum_dev_std = sum_dev_std + ((df.at[i, column_name] - mean_obj) ** 2)
        sigma = mt.sqrt(sum_dev_std / iter_range)

        if objective_function[0] == 'MAX':
            obj = -obj_old
            if limit_exclusions[0] == 'YES' and sigma > dev_std:
                obj = -mean_obj * limit_exclusions[1]


        else:
            obj = obj_old
            if limit_exclusions[0] == 'YES' and sigma > dev_std:
                obj = mean_obj * limit_exclusions[1]
            pass

    except:
        obj = 0



    path_parent = os.path.dirname(os.getcwd())
    os.chdir(path_parent)



    return obj


def pso(func, lb, ub, ieqcons=[], f_ieqcons=None, args=(), kwargs={},
        swarmsize=100, omega=0.5, phip=0.5, phig=0.5, maxiter=100,
        minstep=1e-8, minfunc=1e-8, debug=False, processes=1,
        particle_output=False):
    """
    Perform a particle swarm optimization (PSO)

    Parameters
    ==========
    func : function
        The function to be minimized
    lb : array
        The lower bounds of the design variable(s)
    ub : array
        The upper bounds of the design variable(s)

    Optional
    ========
    ieqcons : list
        A list of functions of length n such that ieqcons[j](x,*args) >= 0.0 in
        a successfully optimized problem (Default: [])
    f_ieqcons : function
        Returns a 1-D array in which each element must be greater or equal
        to 0.0 in a successfully optimized problem. If f_ieqcons is specified,
        ieqcons is ignored (Default: None)
    args : tuple
        Additional arguments passed to objective and constraint functions
        (Default: empty tuple)
    kwargs : dict
        Additional keyword arguments passed to objective and constraint
        functions (Default: empty dict)
    swarmsize : int
        The number of particles in the swarm (Default: 100)
    omega : scalar
        Particle velocity scaling factor (Default: 0.5)
    phip : scalar
        Scaling factor to search away from the particle's best known position
        (Default: 0.5)
    phig : scalar
        Scaling factor to search away from the swarm's best known position
        (Default: 0.5)
    maxiter : int
        The maximum number of iterations for the swarm to search (Default: 100)
    minstep : scalar
        The minimum stepsize of swarm's best position before the search
        terminates (Default: 1e-8)
    minfunc : scalar
        The minimum change of swarm's best objective value before the search
        terminates (Default: 1e-8)
    debug : boolean
        If True, progress statements will be displayed every iteration
        (Default: False)
    processes : int
        The number of processes to use to evaluate objective function and
        constraints (default: 1)
    particle_output : boolean
        Whether to include the best per-particle position and the objective
        values at those.

    Returns
    =======
    g : array
        The swarm's best known position (optimal design)
    f : scalar
        The objective value at ``g``
    p : array
        The best known position per particle
    pf: arrray
        The objective values at each position in p

    """

    assert len(lb) == len(ub), 'Lower- and upper-bounds must be the same length'
    assert hasattr(func, '__call__'), 'Invalid function handle'
    lb = np.array(lb)
    ub = np.array(ub)
    assert [ub[i] > lb[i] for i in range(0, len(ub))], 'All upper-bound values must be greater than lower-bound values'

    vhigh = np.abs(ub - lb)
    vlow = -vhigh

    # Check for constraint function(s) #########################################
    if f_ieqcons is None:
        if not len(ieqcons):
            if debug:
                print('No constraints given.')
            cons = _cons_none_wrapper
        else:
            if debug:
                print('Converting ieqcons to a single constraint function')
            cons = partial(_cons_ieqcons_wrapper, ieqcons, args, kwargs)
    else:
        if debug:
            print('Single constraint function given in f_ieqcons')
        cons = partial(_cons_f_ieqcons_wrapper, f_ieqcons, args, kwargs)
    is_feasible = partial(_is_feasible_wrapper, cons)

    # Initialize the multiprocessing module if necessary
    if processes > 1:
        import multiprocessing
        mp_pool = multiprocessing.Pool(processes)

    # Initialize the particle swarm ############################################
    S = swarmsize
    D = len(lb)  # the number of dimensions each particle has
    x = np.random.rand(S, D)  # particle positions
    # Initialize the particle's position
    x = lb + x * (ub - lb)

    # Imposing input starting point to a particle
    flow, n, objective_function, typ, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel, start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size, varbound = args
    s = 0
    for i in range(0, len(design_variables_alpha)):
        if design_variables_alpha[i] == 'Y':
            x[0][s] = start_alpha[i]
            s = s + 1
    for i in range(0, len(design_variables_dist)):
        for j in range(0, len(design_variables_dist[i])):
            if design_variables_dist[i][j] == 'Y':
                x[0][s] = start_dist[i + 1][j]
                s = s + 1

    for i in range(0, len(design_variables_dist)):
        if design_variables_crel[i] == 'Y':
            x[0][s] = start_crel[i + 1]
            s = s + 1

    k = 0
    for i in range(0, len(typ)):
        if typ[i] == "IGP":
            for j in range(0, len(design_variables_params[k])):
                if design_variables_params[k][j] == 'Y':
                    x[0][s] = start_params[i][j]
                    s = s + 1
        k = k + 1

    v = np.zeros_like(x)  # particle velocities
    p = np.zeros_like(x)  # best particle positions
    fx = np.zeros(S)  # current particle function values
    ID = np.zeros(S)
    fs = np.zeros(S, dtype=bool)  # feasibility of each particle
    fp = np.ones(S) * np.inf  # best particle function values
    g = []  # best swarm position
    fg = np.inf  # best swarm position starting value
    data_id = open('HistoryPSO.txt', 'a')
    if objective_function[0] == 'MIN':
        data_id.write('Minimisation optimisation.\n')
    else:
        data_id.write('Maximisation optimisation.\n')
    data_id.close()
    PROC = []
    # Calculate objective and constraints for each particle
    if processes > 1:
        if processes > len(x) + 1:
            processes = len(x) + 1
            data_id = open('HistoryPSO.txt', 'a')
            data_id.write(
                'Invalid definition for number of processes: maximum number is equal to [(n° of design variables) + 1].\nAs consequence, number of processes is reduced to %d.\n' % processes)
            data_id.close()
        i = 0
        fs = np.array(mp_pool.map(is_feasible, x))
        while i <= len(x):
            for j in range(i, min(processes + i, len(x))):
                pp = mp.Process(target=obj_pso, args=(x, flow, n, objective_function, typ, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel, start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size, varbound, j))
                PROC.append(pp)
                pp.start()
            for proc in PROC:
                proc.join()
            i = i + processes
    else:
        for j in range(0, len(x)):
            obj_pso(x, flow, n, objective_function, typ, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel, start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size, varbound, j)
            fs[j] = is_feasible(x[j, :])

    for i in range(0, len(x)):
            fx[i] = obj_evaluation_PSO(objective_function, limit_exclusions, iter_range, dev_std, i, x[i])
            data_id = open('HistoryPSO.txt', 'a')
            data_id.write('New design variables for process %d: ' % i)
            for j in range(0, len(x[i])):
                if j < len(x[i]) - 1:
                    data_id.write('%.8f, ' % x[i][j])
                else:
                    data_id.write('%.8f.\n' % x[i][j])
            if objective_function[0] == 'MAX':
                data_id.write('Fitness for process %d: %.8f.\n' % (i, -fx[i]))
            else:
                data_id.write('Fitness for process %d: %.8f.\n' % (i, fx[i]))
            data_id.close()

    data_id = open('HistoryPSO.txt', 'a')

    if objective_function[0] == 'MAX':
        data_id.write('Starting configuration fitness: %.8f.\n' % -fx[0])
    else:
        data_id.write('Starting configuration fitness: %.8f.\n' % fx[0])
    data_id.close()

    # Store particle's best position (if constraints are satisfied)
    i_update = np.logical_and((fx < fp), fs)
    p[i_update, :] = x[i_update, :].copy()
    fp[i_update] = fx[i_update]

    # Update swarm's best position
    i_min = np.argmin(fp)
    if fp[i_min] < fg:
        fg = fp[i_min]
        g = p[i_min, :].copy()
        for i in range(0, len(fx)):
            if i != i_min:
                shutil.rmtree('./Process_%d' % i)
            else:
                os.rename(('./Process_%d' % i_min), './Opt_0')
    else:
        # At the start, there may not be any feasible starting point, so just
        # give it a temporary "best" point since it's likely to change
        g = x[0, :].copy()
        for i in range(0, len(fx)):
            shutil.rmtree('./Process_%d' % ID[i])

    # Initialize the particle's velocity
    v = vlow + np.random.rand(S, D) * (vhigh - vlow)

    data_id = open('HistoryPSO.txt', 'a')
    if objective_function[0] == 'MAX':
        data_id.write('Best Fitness (iteration %d): %.8f.\n' % (0, -fg))
    else:
        data_id.write('Best Fitness (iteration %d): %.8f.\n' % (0, fg))
    data_id.write('Best Point (iteration %d): ' % 0)
    for i in range(0, len(g)):
        if i < len(g) - 1:
            data_id.write('%.8f, ' % g[i])
        else:
            data_id.write('%.8f.\n' % g[i])
    data_id.write('Velocity (iteration %d): ' % 0)
    for i in range(0, len(v)):
        if i < len(v) - 1:
            data_id.write('[%.8f, %.8f], ' % (v[i][0], v[i][1]))
        else:
            data_id.write('[%.8f, %.8f].\n' % (v[i][0], v[i][1]))
    data_id.close()

    # Iterate until termination criterion met ##################################
    it = 1
    while it <= maxiter:
        rp = np.random.uniform(size=(S, D))
        rg = np.random.uniform(size=(S, D))
        # Update the particles velocities
        v = omega * v + phip * rp * (p - x) + phig * rg * (g - x)
        # Update the particles' positions
        x = x + v

        # Correct for bound violations
        maskl = x < lb
        masku = x > ub
        x = x * (~np.logical_or(maskl, masku)) + lb * maskl + ub * masku

        # Update objectives and constraints
        if processes > 1:
            if processes > len(x) + 1:
                processes = len(x) + 1
                data_id = open('HistoryPSO.txt', 'a')
                data_id.write(
                    'Invalid definition for number of processes: maximum number is equal to [(n° of design variables) + 1].\nAs consequence, number of processes is reduced to %d.\n' % processes)
                data_id.close()
            i = 0
            fs = np.array(mp_pool.map(is_feasible, x))
            while i <= len(x):
                for j in range(i, min(processes + i, len(x))):
                    pp = mp.Process(target=obj_pso, args=(
                    x, flow, n, objective_function, typ, design_variables_alpha, design_variables_dist,
                    design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel,
                    start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes, y_plus, thick,
                    progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr,
                    semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil,
                    limit_exclusions, dev_std, iter_range, farfield_size, varbound, j))
                    PROC.append(pp)
                    pp.start()
                for proc in PROC:
                    proc.join()
                i = i + processes
        else:
            for j in range(0, len(x)):
                obj_pso(x, flow, n, objective_function, typ, design_variables_alpha, design_variables_dist,
                                design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel,
                                start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes,
                                y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts,
                                wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining,
                                ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size,
                                varbound, j)
                fs[j] = is_feasible(x[j, :])

        for i in range(0, len(x)):
            fx[i] = obj_evaluation_PSO(objective_function, limit_exclusions, iter_range, dev_std, i, x[i])
            data_id = open('HistoryPSO.txt', 'a')
            data_id.write('New design variables for process %d: ' % i)
            for j in range(0, len(x[i])):
                if j < len(x[i]) - 1:
                    data_id.write('%.8f, ' % x[i][j])
                else:
                    data_id.write('%.8f.\n' % x[i][j])
            if objective_function[0] == 'MAX':
                data_id.write('Fitness for process %d: %.8f.\n' % (i, -fx[i]))
            else:
                data_id.write('Fitness for process %d: %.8f.\n' % (i, fx[i]))
            data_id.close()

        # Store particle's best position (if constraints are satisfied)
        i_update = np.logical_and((fx < fp), fs)
        p[i_update, :] = x[i_update, :].copy()
        fp[i_update] = fx[i_update]

        # Compare swarm's best position with global best position
        i_min = np.argmin(fp)
        if fp[i_min] < fg:
            if debug:
                print('New best for swarm at iteration {:}: {:} {:}' \
                      .format(it, p[i_min, :], fp[i_min]))

            p_min = p[i_min, :].copy()
            stepsize = np.sqrt(np.sum((g - p_min) ** 2))

            if np.abs(fg - fp[i_min]) <= minfunc:
                for i in range(0, len(fx)):
                    shutil.rmtree('./Process_%d' % i)
                print('Stopping search: Swarm best objective change less than {:}' \
                      .format(minfunc))
                if particle_output:
                    return p_min, fp[i_min], p, fp
                else:
                    return p_min, fp[i_min]
            elif stepsize <= minstep:
                for i in range(0, len(fx)):
                    shutil.rmtree('./Process_%d' % i)
                print('Stopping search: Swarm best position change less than {:}' \
                      .format(minstep))
                if particle_output:
                    return p_min, fp[i_min], p, fp
                else:
                    return p_min, fp[i_min]


            else:
                g = p_min.copy()
                fg = fp[i_min]
                for i in range(0, len(fx)):
                    if i != i_min:
                        shutil.rmtree('./Process_%d' % i)
                    else:
                        os.rename(('./Process_%d' % i_min), ('./Opt_%d' % it))
                        for j in range(it-1, 0, -1):
                            print(j)
                            try:
                                shutil.rmtree('./Opt_%d' % j)
                            except:
                                pass


        data_id = open('HistoryPSO.txt', 'a')
        if objective_function[0] == 'MAX':
            data_id.write('Best Fitness (iteration %d): %.8f.\n' % (it, -fg))
        else:
            data_id.write('Best Fitness (iteration %d): %.8f.\n' % (it, fg))

        data_id.write('Best Point (iteration %d): ' % it)
        for i in range(0, len(g)):
            if i < len(g) - 1:
                data_id.write('%.8f, ' % g[i])
            else:
                data_id.write('%.8f.\n' % g[i])
        data_id.write('Velocity (iteration %d): ' % it)
        for i in range(0, len(v)):
            if i < len(v) - 1:
                data_id.write('[%.8f, %.8f], ' % (v[i][0], v[i][1]))
            else:
                data_id.write('[%.8f, %.8f].\n' % (v[i][0], v[i][1]))
        data_id.close()

        if debug:
            print('Best after iteration {:}: {:} {:}'.format(it, g, fg))
        it += 1


    print('Stopping search: maximum iterations reached --> {:}'.format(maxiter))
    import matplotlib.pyplot as plt

    if fg >= 0:
        coeff = [0.8, 1.2]
    else:
        coeff = [1.2, 0.8]

    fg_min = coeff[0] * fg
    fg_max = coeff[1] * fg

    plt.figure()
    plt.plot(it, fg, 'k', linewidth=3, label='PSO history')
    plt.ylim(fg_min, fg_max)
    plt.xlim(0, it + 1)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    plt.legend(fontsize=14)
    plt.legend()
    plt.savefig('Optimisation_history.png')
    plt.close()
    if not is_feasible(g):
        print("However, the optimization couldn't find a feasible design. Sorry")
    if particle_output:
        return g, fg, p, fp
    else:
        return g, fg


def obj_pso(x, flow, n, objective_function, typ, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel, start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size, varbound, k):
    path = ('./Process_%d' % k)
    if not os.path.exists(path):
        os.makedirs(path)

    nairfoils = len(typ)
    for i in range(0, len(typ)):
        if typ[i] == ('MY_FILE_%d' % (i + 1)):
            shutil.copy(('MY_FILE_%d.dat' % (i + 1)), path)
    shutil.copy('cfd_config_settings.cfg', path)
    shutil.copy('CFD_run.sh', path)
    shutil.copy('Geo2.sh', path)
    shutil.copy('mesh_gen.py', path)
    shutil.copy('wing_launch_opt.py', path)
    shutil.copy('multiGeometry.py', path)
    shutil.copy('airfoilShape.py', path)
    shutil.copy('camberline.py', path)
    shutil.copy('cfg_printing.py', path)
    shutil.copy('getThickParamIGP.py', path)
    shutil.copy('input_check_opt.py', path)
    shutil.copy('lift_driver_mode.py', path)
    shutil.copy('main_configuration_file.py', path)
    shutil.copy('pts_gen.py', path)
    shutil.copy('wall_spacing.py', path)
    os.chdir(path)
    # with open('readme.txt') as f:
    #     old_N = f.readlines()
    #
    # dummy_id.write(N)
    alpha = [None] * len(start_alpha)
    dist = [None] * len(start_dist)
    crel = [None] * len(start_crel)
    params = [None] * len(start_params)

    dist[0] = [0, 0]
    crel[0] = 1
    s = 0
    for i in range(0, len(design_variables_alpha)):
        if design_variables_alpha[i] == 'Y':
            alpha[i] = x[k][s]
            s = s + 1
        else:
            alpha[i] = start_alpha[i]

    for i in range(0, len(design_variables_dist)):
        dist[i + 1] = [None, None]
        for j in range(0, len(design_variables_dist[i])):
            if design_variables_dist[i][j] == 'Y':
                dist[i + 1][j] = x[k][s]
                s = s + 1
            else:
                dist[i + 1][j] = start_dist[i + 1][j]

    for i in range(0, len(design_variables_crel)):
        if design_variables_crel[i] == 'Y':
            crel[i + 1] = x[k][s]
            s = s + 1
        else:
            crel[i + 1] = start_crel[i + 1]

    for i in range(0, len(design_variables_params)):
        if typ[i] == 'IGP':
            params[i] = [None] * len(design_variables_params[i])
            for j in range(0, len(design_variables_params[i])):
                if design_variables_params[i][j] == 'Y':
                    params[i][j] = x[k][s]
                    s = s + 1
                else:
                    params[i][j] = start_params[i][j]
        else:
            params[i] = [None] * len(start_params[i])
            for j in range(0, len(start_params[i])):
                    params[i][j] = start_params[i][j]


    wing_launch_opt(flow, n, params, alpha, dist, crel, typ, TE, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, farfield_size)

    VAR = [['KIND_TURB_MODEL', KIND_TURB_MODEL], ['KIND_TRANS_MODEL', KIND_TRANS_MODEL],
           ['FREESTREAM_TURBULENCEINTENSITY', FREESTREAM_TURBULENCEINTENSITY], ['RESTART_SOL', RESTART_SOL],
           ['MACH_NUMBER', Mach], ['AOA', alpha[0]], ['SIDESLIP_ANGLE', SIDESLIP_ANGLE],
           ['FREESTREAM_TEMPERATURE', temp], ['GAS_CONSTANT', GAS_CONSTANT],
           ['FLUID_MODEL', FLUID_MODEL], ['GAMMA_VALUE', GAMMA_VALUE], ['REYNOLDS_NUMBER', Re],
           ['REYNOLDS_LENGTH', REYNOLDS_LENGTH], ['VISCOSITY_MODEL', VISCOSITY_MODEL], ['MU_CONSTANT', mu],
           ['MU_T_REF', MU_T_REF], ['SUTHERLAND_CONSTANT', SUTHERLAND_CONSTANT],
           ['CONDUCTIVITY_MODEL', CONDUCTIVITY_MODEL], ['PRANDTL_LAM', PRANDTL_LAM], ['PRANDTL_TURB', PRANDTL_TURB],
           ['REF_ORIGIN_MOMENT_X', REF_ORIGIN_MOMENT_X], ['REF_ORIGIN_MOMENT_Y', REF_ORIGIN_MOMENT_Y],
           ['REF_ORIGIN_MOMENT_Z', REF_ORIGIN_MOMENT_Z], ['REF_LENGTH', REF_LENGTH], ['REF_AREA', REF_AREA],
           ['CFL_NUMBER', CFL_NUMBER], ['CFL_ADAPT', CFL_ADAPT], ['CFL_ADAPT_PARAM', CFL_ADAPT_PARAM],
           ['ITER', ITER],
           ['CONV_NUM_METHOD_FLOW', CONV_NUM_METHOD_FLOW], ['MUSCL_FLOW', MUSCL_FLOW],
           ['SLOPE_LIMITER_FLOW', SLOPE_LIMITER_FLOW], ['VENKAT_LIMITER_COEFF', VENKAT_LIMITER_COEFF],
           ['JST_SENSOR_COEFF', JST_SENSOR_COEFF], ['CONV_FIELD', CONV_FIELD],
           ['CONV_RESIDUAL_MINVAL', CONV_RESIDUAL_MINVAL], ['CONV_STARTITER', CONV_STARTITER],
           ['CONV_CAUCHY_ELEMS', CONV_CAUCHY_ELEMS], ['CONV_CAUCHY_EPS', CONV_CAUCHY_EPS],
           ['MESH_FILENAME', 'Test.su2'], ['SCREEN_OUTPUT', SCREEN_OUTPUT], ['HISTORY_OUTPUT', HISTORY_OUTPUT],
           ['OUTPUT_FILES', OUTPUT_FILES], ['TABULAR_FORMAT', TABULAR_FORMAT],
           ['SCREEN_WRT_FREQ_INNER', SCREEN_WRT_FREQ_INNER], ['OUTPUT_WRT_FREQ', OUTPUT_WRT_FREQ]]

    VAR2 = [['FIXED_CL_MODE', FIXED_CL_MODE], ['TARGET_CL', TARGET_CL], ['DCL_DALPHA', DCL_DALPHA],
            ['UPDATE_AOA_ITER_LIMIT', UPDATE_AOA_ITER_LIMIT], ['ITER_DCL_DALPHA', ITER_DCL_DALPHA],
            ['EVAL_DOF_DCX', EVAL_DOF_DCX]]

    cfg_printing(nairfoils, VAR, VAR2, flow) #, k)
    # if not os.path.exists('./Candidate_%d' % i):
    #     os.makedirs('./Candidate_%d' % i)
    subp = subprocess.call("./CFD_run.sh")

    path_parent = os.path.dirname(os.getcwd())
    os.chdir(path_parent)

    return []


def opt_method_pso(flow, n, objective_function, typ, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel, start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size, varbound):

    swarmsize, omega, phip, phig, maxiter, minstep, minfunc, debug, processes, particle_output = main_PSO_algo_params()
    if swarmsize == 'DEFAULT':
        swarmsize = 20
    if omega == 'DEFAULT':
        omega = 0.5
    if phip == 'DEFAULT':
        phip = 0.5
    if phig == 'DEFAULT':
        phig = 0.5
    if maxiter == 'DEFAULT':
        maxiter = 100
    if minstep == 'DEFAULT':
        minstep = 1e-5
    if minfunc == 'DEFAULT':
        minfunc = 1e-3
    if debug == 'DEFAULT':
        debug = False
    if processes == 'DEFAULT':
        processes = 1
    if particle_output == 'DEFAULT':
        particle_output = False

    data_id = open('HistoryPSO.txt', 'w+')
    data_id.write('-- PSO OPTIMISATION -- \n\n')
    data_id.write(
        'swarmsize: %d.\nomega: %.5f.\nphip: %.5f.\nphig: %.5f.\nmaxiter: %d.\nminstep: %.10f.\nminfunc: %.10f.\nprocesses: %d\n.debug: %s.\nparticle_output: %s.\n' % (
        swarmsize, omega, phip, phig, maxiter, minstep, minfunc, processes, str(debug), str(particle_output)))

    lb = []
    ub = []

    for i in range(0, len(varbound)):
        lb.append(varbound[i][0])
        ub.append(varbound[i][1])


    data_id.write('Lower bounds: ')
    for i in range(0, len(lb)):
        if i < len(lb) - 1:
            data_id.write('%.5f, ' % lb[i])
        else:
            data_id.write('%.5f.\n' % lb[i])

    data_id.write('Upper bounds: ')
    for i in range(0, len(ub)):
        if i < len(ub) - 1:
            data_id.write('%.5f, ' % ub[i])
        else:
            data_id.write('%.5f.\n' % ub[i])
    data_id.close()
    lb = (np.asarray(lb)).tolist()
    ub = (np.asarray(ub)).tolist()

    args = (flow, n, objective_function, typ, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, TE, start_dist, start_alpha, start_crel, start_params, bounds_alpha, bounds_dist, bounds_crel, bounds_params, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, limit_exclusions, dev_std, iter_range, farfield_size, varbound)

    xopt, fopt = pso(obj_pso, lb, ub, args=args, swarmsize=swarmsize, omega=omega, phip=phip, phig=phig,
                     maxiter=maxiter, minstep=minstep, minfunc=minfunc, debug=debug, processes=processes,
                     particle_output=particle_output)

    print('The optimum is at:')
    print('    {}'.format(xopt))
    print('Optimal function values:')
    print('    Optimal Fitness         : {}'.format(fopt))


    return []
