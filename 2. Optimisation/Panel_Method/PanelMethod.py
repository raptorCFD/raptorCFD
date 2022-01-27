from Panel_Method.resultscheck import resultscheck
import time
from Panel_Method.VelocityField import velocityField
from airfoilShape import *
from Panel_Method.solverHS import solverHS

def plotcandidate(nairfoils, x, y, index, uField, vField, Xgrid, Ygrid, k, showNodes, typ, Cl, Cd):
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    w = 3
    speed = np.sqrt(uField ** 2 + vField ** 2)

    fig = plt.figure(figsize=(7, 9))
    gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1, 1, 2])
    ax0 = fig.add_subplot(gs[0, 0])
    ax0.streamplot(Xgrid, Ygrid, uField, vField, density=[0.5, 1])

    plt.ylim((-1.5, 1.5))
    strCl = []
    strCd = []

    for j in range(0, nairfoils):
        if typ[j] == ('MY_FILE_%d' % (j+1)):
            x_plot = [0] * index[j]
            y_plot = [0] * index[j]
            for p in range(0, index[j]):
                x_plot[k] = x[j, p]
                y_plot[k] = y[j, p]
            plt.plot(x_plot, y_plot, '-', color=np.random.rand(3,1), linewidth=2)
            if showNodes == 1:
                plt.scatter(x_plot, y_plot, '.', 'black')
        else:
            plt.plot(x, y, '-', color=np.random.rand(3,1), linewidth=2)
            if showNodes == 1:
                plt.scatter(x_plot, y_plot, '.', 'black')

            strCl.append(str(Cl(k, j)))
            strCd.append(str(Cd(k, j)))

            if j != nairfoils:
                strCl.append('   ')
                strCd.append('   ')

            plt.title('\nCL (left to right):   ' + strCl + '\nCD (left to right):   ' + strCd + '\n\n')
    plt.axis('equal')
    plt.savefig('HS_PreliminaryCandidate_%d.png' % k)
    plt.close()

    return []


def HS_panelmethod(n, start_alpha, start_dist, start_crel, start_params, design, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, ref_airfoil, typ, TE, obj_function, bounds_margin, varbound):

    from main_preliminary_study import HS
    npoint, showNodes, obj = HS()
    npoint = (0.5 * npoint) + 1     # to obtain the desired input, due to multiGeometry.py structure
    nairfoils = len(start_alpha)
    N1 = len(design)
    N2 = len(design[0])
    alpha, dist, crel, params = [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [None] * nairfoils
    xmax, ymax, Cl, Cd, maxdCp, f, TOT_Cl, TOT_Cd = [None] * N1, [None] * N1, [None] * N1, [None] * N1, [None] * N1, [None] * N1, [None] * N1, [None] * N1

    dist[0] = [0, 0]
    crel[0] = 1

    for k in range(0, N1):
        start = time.process_time()

        print('Hess-Smith panel method preliminary study: candidate %d analysis...' % k)
        s = 0
        for i in range(0, len(design_variables_alpha)):
            if design_variables_alpha[i] == 'Y':
                alpha[i] = design[k][s]
                s = s + 1
            else:
                alpha[i] = start_alpha[i]

        for i in range(0, len(design_variables_dist)):
            dist[i + 1] = [None, None]
            for j in range(0, len(design_variables_dist[i])):
                if design_variables_dist[i][j] == 'Y':
                    dist[i + 1][j] = design[k][s]
                    s = s + 1
                else:
                    dist[i + 1][j] = start_dist[i + 1][j]

        for i in range(0, len(design_variables_crel)):
            if design_variables_crel[i] == 'Y':
                crel[i + 1] = design[k][s]
                s = s + 1
            else:
                crel[i + 1] = start_crel[i + 1]

        for i in range(0, len(design_variables_params)):
            if typ[i] == "IGP":
                reference = len(design_variables_params[i])
            else:
                reference = len(start_params[i])

            params[i] = [None] * reference

            for j in range(0, reference):
                if design_variables_params[i][j] == 'Y':
                    params[i][j] = design[k][s]
                    s = s + 1
                else:
                    params[i][j] = start_params[i][j]

        # Now the geometry is defined. We can proceed to the specific HS panel method calculus for candidate "design[k]"
        Cl[k], Cd[k], xmax[k], ymax[k], maxdCp[k], x, y, index, p, p1, SOL, metaPan, nairfoils = solverHS(alpha, npoint, params, ref_airfoil, typ, TE, dist, crel, n)

        #uField, vField, Xgrid, Ygrid = velocityField(p, nairfoils, alpha, 1, SOL, x, y, metaPan)

        TOT_Cl[k] = sum(Cl[k])
        TOT_Cd[k] = sum(Cd[k])

        # Plot candidate n° k
        #plotcandidate(nairfoils, x, y, index, uField, vField, Xgrid, Ygrid, k, showNodes, typ, Cl, Cd)

        print('Candidate %d requested time for Hess-Smith panel method evaluation: (%.3f s).\n' % (k, time.process_time() - start))


    # A final check about computed results is executed, trying to also reduce the design space by reducing lower and upper bounds.
    alpha_opt_preliminary, dist_opt_preliminary, crel_opt_preliminary, params_opt_preliminary, design_opt, ub, lb, OBJ, bounds = resultscheck(design, obj_function, obj, maxdCp, TOT_Cl, TOT_Cd, start_alpha, start_dist, start_crel, start_params, design_variables_alpha, design_variables_dist, design_variables_crel, design_variables_params, N1, N2, bounds_margin, varbound)

    return alpha_opt_preliminary, crel_opt_preliminary, dist_opt_preliminary, params_opt_preliminary, ub, lb, OBJ, bounds
