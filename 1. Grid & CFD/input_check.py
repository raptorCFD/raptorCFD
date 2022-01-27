import sys


def input_check(params, alpha, dist, crel, typ, TE, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, flow, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, n, farfield_size):

    print('Checking consistency of given input values...')


    # Chosen parametrization
    for i in range(0, len(typ)):
        if typ[i] != 'IGP' and typ[i] != 'NACA4' and typ[i] != 'NACA5' and typ[i] != 'NACA4_MOD' and typ[i] != 'NACA16' and typ[i] != 'BENZING' and typ[i] != ('MY_FILE_%d' % (i+1)):
            print(
                "Error: the chosen parametrization should be one of the following: IGP, NACA4, NACA5, NACA4_MOD, NACA16 or MY_FILE_%d for airfoil n° %d. \n The written should be as shown above. Check the item: 'typ'." % (i+1, i+1))
            sys.exit()
            # Option "BICONVEX" to be developed. "BENZING" to be completed

    # Length of each parameter
    if len(params) != len(typ):
        print(
            "Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \nCheck the following item length: 'params'.")
        sys.exit()
    elif len(alpha) != len(typ):
        print(
            "Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \n Check the following item length: 'alpha'.")
        sys.exit()
    elif len(dist) != len(typ):
        print(
            "Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \n Check the following item length: 'dist'.")
        sys.exit()
    elif len(TE) != len(typ):
        print(
            "Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \n Check the following item length: 'TE'.")
        sys.exit()
    elif len(crel) != len(typ):
        print(
            "Error: length of each item should be equal to the length of item (" + str(len(typ)) + "). \n Check the following item length: 'crel'.")
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

    # Parameters input check
    for i in range(0, len(params)):
        for j in range(0, len(params[i])):
            if params[i][j] != ("MY_FILE_%d" % (i+1)) and type(params[i][j]) != int and type(params[i][j]) != float:
                print("Error: each value inserted in 'params' must be an integer or the written 'MY_FILE_%d' for airfoil n° %d. \n Check the item: 'params'." % (i+1, i+1))
                sys.exit()
            else:
                pass

    # Distance between airfoils item
    for i in range(0, len(typ)):
        if type(alpha[i]) != int and type(alpha[i]) != float:
            print(
                "Error: angle of attack of each airfoil should be written as integer or float number. \n Check the item: 'alpha'.")
            sys.exit()
        elif (type(alpha[i]) == int or type(alpha[i]) == float) and (alpha[i] > 90 or alpha[i] < -90):
            print(
                "Error: angle of attack of each airfoil must be ranging from -90° to +90°. \n Check the item: 'alpha'.")
            sys.exit()


        if type(crel[i]) != int and type(crel[i]) != float:
            print("Error: relative chord of each airfoil should be defined by an integer or float number greater than 0. \n Check the item: 'crel'.")
            sys.exit()
        elif (type(crel[i]) != int or type(crel[i]) != float) and crel[i] <= 0:
            print("Error: relative chord of each airfoil should be defined by an integer or float number greater than 0. \n Check the item: 'crel'.")
            sys.exit()

        for j in range(0, len(dist[i])):
            if type(dist[i][j]) != int and type(dist[i][j]) != float:
                print(
                    "Error: distance between airfoils should be an integer or float number. \n Check the item: 'dist'.")
                sys.exit()
            else:
                pass

    for k in range(0, len(typ)):
        if TE[k] != ('MY_FILE_%d' % (k+1)) and TE[k] != 'Y' and TE[k] != 'N':
            print("Error: each term inserted in 'TE' item must be 'Y', 'N' or 'MY_FILE_%d' for airfoil n° %d. \n Check the item: 'TE'." % (k+1, k+1))
            sys.exit()
        else:
            pass

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
    if (type(nodes) == int or type(nodes) == float) and nodes <= 0:
        print("Error: the number of airfoils must be an integer or a float number greater than 0.")
        sys.exit()
    elif type(nodes) != int and type(nodes) != float:
        print("Error: the number of airfoils must be an integer or a float number greater than 0.")
        sys.exit()
    else:
        pass

    for i in range(0, len(y_plus)):
        if type(y_plus[i]) != int and type(y_plus[i]) != float:
            print('Error: dimensionless wall distance "y_plus" must be a vector of numerical values, greater than 0 and lower than 1. \n Check item "y_plus".')
            sys.exit()
        elif (type(y_plus[i]) == int or type(y_plus[i]) == float) and (y_plus[i] <= 0 or y_plus[i] > 1):
            print('Error: dimensionless wall distance "y_plus" must be a vector of numerical values, greater than 0 and lower than 1. \n Check item "y_plus".')
            sys.exit()
            pass


    for i in range(0, len(thick)):
        if (type(thick[i]) != int and type(thick[i]) != float):
            print('Error: structured region thickness must be a vector of numerical values, each greater than 0. Pay attention to avoid intersections between structured regions and airfoils\' geometries. \n Check item "thick".')
            sys.exit()
        elif (type(thick[i]) == int or type(thick[i]) == float) and thick[i] <= 0:
            print('Error: structured region thickness must be a vector of numerical values, each greater than 0. Pay attention to avoid intersections between structured regions and airfoils\' geometries. \n Check item "thick".')
            sys.exit()

    for i in range(0, len(progr)):
        if ((type(progr[i]) != int and type(progr[i]) != float)):
            print('Error: structured region progression item must be a vector of numerical values, each one greater than 1 and lower than 1.3. \n Check item "progr".')
            sys.exit()
        elif ((type(progr[i]) == int or type(progr[i]) == float)) and (progr[i] < 1 or progr[i] > 1.3):
            print('Error: structured region progression item must be a vector of numerical values, each one greater than 1 and lower than 1.3. \n Check item "progr".')
            sys.exit()


    if (type(temp) != int and type(temp) != float):
        print('Error: free-stream temperature item must be a numerical value greater than 0. \n Check item "temp".')
        sys.exit()
    elif (type(temp) == int or type(temp) == float) and temp <= 0:
        print('Error: free-stream temperature item must be a numerical value greater than 0. \n Check item "temp".')
        sys.exit()

    if (type(Mach) != int and type(Mach) != float):
        print('Error: free-stream Mach item must be a numerical value greater than 0. \n Check item "Mach".')
        sys.exit()
    elif (type(Mach) == int or type(Mach) == float) and Mach <= 0:
        print('Error: free-stream Mach item must be a numerical value greater than 0. \n Check item "Mach".')
        sys.exit()


    if type(ellipse_dimension) != float or (type(ellipse_dimension) == float and ellipse_dimension < 0.1 and ellipse_dimension > 0.95):
        print('Error: dimension of the ellipse should be defined by a float number between 0.1 and 0.95. \n Check item "ellipse_dimension.')
        sys.exit()

    if (type(Mesh_Algo) != int and Mesh_Algo != 'DEFAULT') or (type(Mesh_Algo) == int and (Mesh_Algo <= 0 or Mesh_Algo >= 10)):
        print('Error: mesh algorithm item accepts only integer values between 1 and 9, or the string "DEFAULT" which corresponds to integer 6 (Frontal-Delaunay). \n Check item "Mesh_Algo".')
        sys.exit()

    if external_pts != 'YES' and external_pts != 'NO':
        print('Error: ellipse\'s external points\' item "external_pts" accepts only strings "YES" or "NO". \n Check item "external_pts".')
        sys.exit()

    if type(wake_length) != int and type(wake_length) != float:
        print('Error: wake length item should be an numerical value greater than 1. It is also suggested a value lower than 250-500 approximately (uncertain since it strictly depends on dimension of geometry chosen). \n Check item "wake_length".')
        sys.exit()
    elif (type(wake_length) == int or type(wake_length) == float) and wake_length < 1:
        print('Error: wake length item should be an numerical value greater than 1. It is also suggested a value lower than 250-500 approximately (uncertain since it strictly depends on dimension of geometry chosen). \n Check item "wake_length".')
        sys.exit()

    if type(wake_progr) != int and type(wake_progr) != float:
        print('Error: wake element progression item should be an numerical value greater than or equal to 1. \n Check item "wake_progr".')
        sys.exit()
    elif (type(wake_progr) == int or type(wake_progr) == float) and wake_progr < 1:
        print('Error: wake element progression item should be an numerical value greater than or equal to 1. \n Check item "wake_progr".')
        sys.exit()

    if type(semicircle_dimension) != int and type(semicircle_dimension) != float:
        print('Error: external semi-circle item should be an numerical value greater than 1. \n Check item "semicircle_dimension".')
        sys.exit()
    elif (type(semicircle_dimension) == int or type(semicircle_dimension) == float) and semicircle_dimension <= 1:
        print('Error: external semi-circle item should be an numerical value greater than 1. \n Check item "semicircle_dimension".')
        sys.exit()

    if type(semicircle_elem_factor) != int and type(semicircle_elem_factor) != float:
        print('Error: external semi-circle item should be an numerical value greater than 1. \n Check item "semicircle_elem_factor".')
        sys.exit()
    elif (type(semicircle_elem_factor) == int or type(semicircle_elem_factor) == float) and semicircle_elem_factor <= 1:
        print('Error: external semi-circle item should be an numerical value greater than 1. \n Check item "semicircle_elem_factor".')
        sys.exit()

    for i in range(0, len(wall_refining)):
        if (type(wall_refining[i]) != int and type(wall_refining[i]) != float):
            print('Error: wall refining\'s item should be an numerical value greater than 0. \n Check item "wall_refining".')
            sys.exit()
        elif (type(wall_refining[i]) == int or type(wall_refining[i]) == float) and wall_refining[i] <= 0:
            print('Error: wall refining\'s item should be an numerical value greater than 0. \n Check item "wall_refining".')
            sys.exit()

    if (type(ellipse_refining) != int and type(ellipse_refining) != float) or ellipse_refining <= 0:
        print('Error: ellipse refining\'s item should be an numerical value greater than 0. \n Check item "ellipse_refining".')
        sys.exit()
    elif (type(ellipse_refining) == int or type(ellipse_refining) != float) and ellipse_refining <= 0:
        print('Error: ellipse refining\'s item should be an numerical value greater than 0. \n Check item "ellipse_refining".')
        sys.exit()

    if (type(ref_airfoil) != int and ref_airfoil != "DEFAULT"):
        print(
            'Error: reference airfoil\'s item should be an integer value greater or equal to 1, or string written "DEFAULT". \n Check item "ref_airfoil".')
        sys.exit()
    elif type(ref_airfoil) == int and (ref_airfoil < 1 or ref_airfoil > len(typ)):
        print(
            'Error: reference airfoil\'s item should be an integer value greater or equal to 1, or string written "DEFAULT". \n Check item "ref_airfoil".')
        sys.exit()

    if (type(n) != int and n != "DEFAULT"):
        print('Error: number of generating points for each airfoil should be an integer value greater than 0, or string written "DEFAULT". \n Check item "n".')
        sys.exit()
    elif type(n) == int and n <= 0:
        print('Error: number of generating points for each airfoil should be an integer value greater than 0, or string written "DEFAULT". \n Check item "n".')
        sys.exit()

    if flow != 'EULER' and flow != 'VISCOUS':
        print('Error: flow type item "flow" should be string written "EULER" or "VISCOUS". \n Check item "flow".')
        sys.exit()

    if type(farfield_size) != int or (type(farfield_size) == int and farfield_size <= sum(crel)):
        print('Error: item "farfield_size" must be an integer value higher than sum of relative chords. \n Check item "farfield_size".')
        sys.exit()

    print('Input check completed! \n')


    return []
