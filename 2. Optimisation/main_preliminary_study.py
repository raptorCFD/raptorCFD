def BB():
    # ================= Box-Behnken (BB)
    # centre: number of points to include. If user inserts "DEFAULT", then a pre-determined number of points are automatically included.
    centre = 'DEFAULT'
    size = 0.025
    processes = 3
    bounds_margin = 20 

    return centre, size, processes, bounds_margin


# def CC():
#     # ================= Central Composite
#     # centre: 2-tuple of center points(one for the factorial block, one for the star block).
#     # alpha: either “orthogonal” (or “o”) or “rotatable” (or “r”).
#     # face: either “circumscribed” (or “ccc”), “inscribed” (or “cci”), or “faced” (or “ccf”).
#     centre = (0, 1)
#     alpha = 'o'
#     face = 'cci'
#     size = 0.025
#
#     return centre, alpha, face, size


def LHS():
    # ================= Latin-hypercube designs
    # samples: an integer that designates the number of sample points to generate for each factor ('DEFAULT' means equal to design variables)
    # criterion: a string that tells lhs how to sample the points.
    #   'None': simply randomizes the points within the intervals.
    #   'center' or 'c': center the points within the sampling intervals
    #   'maximin' or 'm': maximize the minimum distance between points, but place the point in a randomized location within its interval
    #   'centermaximin' or 'cm': same as “maximin”, but centered within the intervals
    #   'correlation' or 'corr': minimize the maximum correlation coefficient
    samples = 10
    criterion = 'c'
    stdvs = [2.5, 3.5] #, 0.0075, 0.002, 0.075, 0.15, 0.09, 0.075, 0.035, 0.05, 0.05, 0.1, 0.175, 0.15, 0.09, 0.075, 0.035, 0.05, 0.05, 0.1, 0.175]  # In sequence 'Y' design variables inserted in main_opt with the following order: alpha, dist, crel, params
    processes = 10
    bounds_margin = 3.1  # Min/Max value of OBJ - VALAREZO-CHIN

    return samples, criterion, stdvs, processes, bounds_margin


def HS():
    # Hess-Smith Panel Method inputs
    npoint = 80
    showNodes = 0  # 1: yes, 0: no
    obj = 'OBJ'  # 'VALAREZO-CHIN' option is valid only for max lift optimisation problems.

    return npoint, showNodes, obj
