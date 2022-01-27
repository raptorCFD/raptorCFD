def resultscheck(design, obj_function, obj, maxdCp, TOT_Cl, TOT_Cd, start_alpha, start_dist, start_crel, start_params,
                 design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, N1, N2,
                 bounds_margin, varbound):
    import numpy as np
    if obj_function[1] == "LIFT":

        if obj == "VALAREZO-CHIN":
            sum_maxdCp = [None] * len(maxdCp)
            for i in range(0, len(maxdCp)):
                if np.any(maxdCp[i]) < 13.5:
                    sum_maxdCp[i] = sum(maxdCp[i])
                else:
                    sum_maxdCp[i] = "VALAREZO-CHIN LIMIT OVERCOME"
                    print('Valarezo-Chin limit on a single element (maxdCp = 13.5) overcome.\n')
            obj_opt = max(sum_maxdCp)
            design_opt = np.argmax(sum_maxdCp)
        elif obj_function[1] == "MAX":
            obj_opt = max(TOT_Cl)
            design_opt = np.argmax(TOT_Cl)
        else:
            obj_opt = min(TOT_Cl)
            design_opt = np.argmax(TOT_Cl)

    elif obj_function[1] == "DRAG":
        obj_opt = min(TOT_Cd)
        design_opt = np.argmin(TOT_Cd)
    else:
        TOT_Eff = [0] * len(TOT_Cl)
        for i in range(0, len(TOT_Cl)):
            TOT_Eff[i] = TOT_Cl[i] / TOT_Cd[i]

            if obj_function[1] == "MAX":
                obj_opt = max(TOT_Eff)
                design_opt = np.argmax(TOT_Eff)
            else:
                obj_opt = min(TOT_Eff)
                design_opt = np.argmin(TOT_Eff)

    print('Optimal preliminary configuration is number %d with fitness %.5f.\n' % (design_opt, obj_opt))
    alpha_opt_preliminary = start_alpha
    dist_opt_preliminary = start_dist
    crel_opt_preliminary = start_crel
    params_opt_preliminary = start_params

    s = 0
    for j in range(0, len(alpha_opt_preliminary)):
        if design_variables_alpha[j] == 'Y':
            alpha_opt_preliminary[j] = design[design_opt][s]
            s = s + 1
        else:
            alpha_opt_preliminary[j] = start_alpha[j]

    for j in range(1, len(dist_opt_preliminary)):
        for k in range(0, len(dist_opt_preliminary[j])):
            if design_variables_dist[j - 1][k] == 'Y':
                dist_opt_preliminary[j][k] = design[design_opt][s]
                s = s + 1
            else:
                dist_opt_preliminary[j][k] = start_dist[j][k]

    for j in range(1, len(crel_opt_preliminary)):
        if design_variables_crel[j-1] == 'Y':
            crel_opt_preliminary[j] = design[design_opt][s]
            s = s + 1
        else:
            crel_opt_preliminary[j] = start_crel[j]

    for j in range(0, len(params_opt_preliminary)):
        for k in range(0, len(params_opt_preliminary[j])):
            if design_variables_params[j - 1][k] == 'Y':
                params_opt_preliminary[j][k] = design[design_opt][s]
                s = s + 1
            else:
                params_opt_preliminary[j][k] = start_params[j][k]

    ub = [-1000] * N2
    lb = [1000] * N2

    if obj_function[1] == "LIFT":
        if obj == "VALAREZO-CHIN":
            REF = sum_maxdCp
        else:
            REF = TOT_Cl

    elif obj_function[1] == "DRAG":
        REF = TOT_Cd
    elif obj_function[1] == "EFFICIENCY":
        REF = [None] * len(TOT_Cl)
        for ii in range(0, len(TOT_Cl)):
            REF[ii] = TOT_Cl[ii] / TOT_Cd[ii]

    if obj_function[1] == "LIFT" and obj == "VALAREZO-CHIN" and obj_opt < bounds_margin:
        print(
            'The optimal fitness found by preliminary Panel Method study is lower then "bounds_margin" given in input.\n')
        print(
            'As consequence, the optimal fitness will be taken as reference, and bounds dictated by samples respective fitness higher than half of optimum.\n')
        for j in range(0, N2):
            for i in range(0, N1):
                if REF[i] > obj_opt * 0.5:
                    if design[i][j] > ub[j]:
                        ub[j] = design[i][j]
                    if design[i][j] < lb[j]:
                        lb[j] = design[i][j]
    else:
        limiter = bounds_margin
        for j in range(0, N2):
            for i in range(0, N1):
                if REF[i] > limiter:
                    if design[i][j] > ub[j]:
                        ub[j] = design[i][j]
                    if design[i][j] < lb[j]:
                        lb[j] = design[i][j]
    for i in range(0, N2):
        if lb[i] == design[design_opt][i]:
            if lb[i] > 0:
                lb[i] = design[design_opt][i] * 0.5
            elif lb[i] < 0:
                lb[i] = design[design_opt][i] * 1.5
            else:
                lb[i] = lb[i] - abs(ub[i]) * 0.5
        if ub[i] == design[design_opt][i]:
            if ub[i] > 0:
                ub[i] = design[design_opt][i] * 1.5
            elif ub[i] < 0:
                ub[i] = design[design_opt][i] * 0.5
        else:
            ub[i] = ub[i] + abs(lb[i]) * 0.5

    # Bounds' check
    bounds = [0] * N2
    for i in range(0, N2):
        bounds[i] = [0, 0]
        bounds[i][0] = lb[i]
        bounds[i][1] = ub[i]
        if bounds[i][0] < varbound[i][0]:
            bounds[i][0] = varbound[i][0]
        if bounds[i][1] > varbound[i][1]:
            bounds[i][1] = varbound[i][1]

    return alpha_opt_preliminary, dist_opt_preliminary, crel_opt_preliminary, params_opt_preliminary, design_opt, ub, lb, REF, bounds
