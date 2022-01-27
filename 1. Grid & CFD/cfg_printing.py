def cfg_printing(nairfoils, var, var2, flow):

    print('Modifying Configuration file...\n')

    fout = open('cfd_config_settings.cfg', "wt")
    fout.write('% Physical governing equations (EULER, NAVIER_STOKES, \n% WAVE_EQUATION, HEAT_EQUATION, FEM_ELASTICITY, \n% POISSON_EQUATION)\n')
    if flow == 'EULER':
        pass
    else:
        flow = 'RANS'
    fout.write('SOLVER= %s\n\n' % flow)

    for i in range(0, len(var)):
        fout.write('%s = %s\n\n' % (var[i][0], var[i][1]))

    fout.write('% Mathematical problem(DIRECT, CONTINUOUS_ADJOINT)\n')
    fout.write('MATH_PROBLEM = DIRECT\n')
    fout.write('% Init option to choose between Reynolds(default) or thermodynamics quantities\n')
    fout.write('% for initializing the solution (REYNOLDS, TD_CONDITIONS)\n')
    fout.write('INIT_OPTION = REYNOLDS\n')
    fout.write('% Free - stream option to choose between density and temperature(default) for\n')
    fout.write('% initializing the solution(TEMPERATURE_FS, DENSITY_FS)\n')
    fout.write('FREESTREAM_OPTION = TEMPERATURE_FS\n')
    fout.write('% Compressible flow non - dimensionalization(DIMENSIONAL, FREESTREAM_PRESS_EQ_ONE,\n')
    fout.write('% FREESTREAM_VEL_EQ_MACH, FREESTREAM_VEL_EQ_ONE)\n')
    fout.write('REF_DIMENSIONALIZATION = DIMENSIONAL\n')

    fout.write('% --------------------------- LIFT DRIVER MODE --------------------------------%\n')
    for i2 in range(0, len(var2)):
        fout.write('%s = %s\n\n' % (var2[i2][0], var2[i2][1]))

    fout.write('\n% THE FOLLOWING VARIABLES WERE NOT DEFINED IN INPUT BY USER, BUT DEFINED AS DEFAULT.\n')
    fout.write('% -------------------- TURBULENT NUMERICAL METHOD DEFINITION ------------------%\n')
    fout.write('% Convective numerical method (SCALAR_UPWIND)\n')
    fout.write('CONV_NUM_METHOD_TURB = SCALAR_UPWIND\n')
    #
    fout.write('% ------------- COMMON PARAMETERS DEFINING THE NUMERICAL METHOD ---------------%\n')
    fout.write('% Numerical method for spatial gradients (GREEN_GAUSS, WEIGHTED_LEAST_SQUARES)\n')
    fout.write('NUM_METHOD_GRAD = WEIGHTED_LEAST_SQUARES\n')
    fout.write('% Time discretization (RUNGE-KUTTA_EXPLICIT, EULER_IMPLICIT, EULER_EXPLICIT)\n')
    fout.write('TIME_DISCRE_FLOW = EULER_IMPLICIT\n')

    fout.write('% ------------------------ LINEAR SOLVER DEFINITION ---------------------------%\n')
    fout.write('% Linear solver for the implicit (or discrete adjoint) formulation (BCGSTAB, FGMRES)\n')
    fout.write('LINEAR_SOLVER= FGMRES\n')
    fout.write('% Preconditioner of the Krylov linear solver (NONE, JACOBI, LINELET)\n')
    fout.write('LINEAR_SOLVER_PREC= ILU\n')
    fout.write('% Min error of the linear solver for the implicit formulation\n')
    fout.write('LINEAR_SOLVER_ERROR= 1E-8\n')
    fout.write('% Max number of iterations of the linear solver for the implicit formulation\n')
    fout.write('LINEAR_SOLVER_ITER= 10\n')

    fout.write('% -------------------- BOUNDARY CONDITION DEFINITION --------------------------%\n')

    if flow == 'RANS':
        fout.write('% Navier-Stokes wall boundary marker(s) (NONE = no marker)\n')
        fout.write('MARKER_HEATFLUX= ( airfoil1, 0.0')
        for j in range(1, nairfoils):
            if j == nairfoils - 2:
                fout.write(', airfoil%d, 0.0' % (j + 1))
            else:
                fout.write(', airfoil%d, 0.0 )\n\n' % (j + 1))
    else:
        fout.write('MARKER_EULER= ( airfoil1, 0.0')
        for j in range(1, nairfoils):
            if j == nairfoils - 2:
                fout.write(', airfoil%d, 0.0' % (j + 1))
            else:
                fout.write(', airfoil%d, 0.0 )\n\n' % (j + 1))


    fout.write('% Marker(s) of the surface to be plotted or designed\n')
    fout.write('MARKER_PLOTTING= ( airfoil1')
    for j in range(1, nairfoils):
        if j == nairfoils - 2:
            fout.write(', airfoil%d,' % (j + 1))
        else:
            fout.write(', airfoil%d )\n\n' % (j + 1))

    fout.write('% Marker(s) of the surface where the functional (Cd, Cl, etc.) will be evaluated\n')
    fout.write('MARKER_MONITORING= ( airfoil1')
    for j in range(1, nairfoils):
        if j == nairfoils - 2:
            fout.write(', airfoil%d,' % (j + 1))
        else:
            fout.write(', airfoil%d )\n\n' % (j + 1))

    fout.write('% Farfield boundary marker(s) (NONE = no marker)\n')
    fout.write('MARKER_FAR= ( farfield )\n\n')

    fout.write('% ------------------------- INPUT / OUTPUT INFORMATION --------------------------- %\n')
    fout.write('% Mesh input file format(SU2, CGNS, NETCDF_ASCII)\n')
    fout.write('MESH_FORMAT = SU2\n\n')

    fout.write('% Mesh output file\n')
    fout.write('MESH_OUT_FILENAME = mesh_out.su2\n\n')

    fout.write('% Restart flow input file\n')
    fout.write('SOLUTION_FILENAME = solution_flow.dat\n\n')

    fout.write('% Output file convergence history(w / o extension)\n')
    fout.write('CONV_FILENAME = history\n\n')

    fout.write('% Output file restart flow\n')
    fout.write('RESTART_FILENAME = restart_flow.dat\n\n')

    fout.write('% Output file flow(w / o extension) variables\n')
    fout.write('VOLUME_FILENAME = flow\n\n')

    fout.write('% Output file surface flow coefficient(w / o extension)\n')
    fout.write('SURFACE_FILENAME = surface_flow\n\n')

    fout.write('WRT_FORCES_BREAKDOWN = YES\n')

    fout.write('% Output file with the forces breakdown\n')
    fout.write('BREAKDOWN_FILENAME = forces_breakdown.dat\n\n')

    fout.write('% Writing frequency for history output\n')
    fout.write('HISTORY_WRT_FREQ_INNER = 1\n\n')


    fout.close()


    print('Configuration file setted! \n') # Check .cfg file for correct settings.\n')

    return []
