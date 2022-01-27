import numpy as np
from airfoilShape import airfoilIGP
from airfoilShape import airfoilNACA4
from airfoilShape import airfoilNACA5
from airfoilShape import airfoilNACA4_MOD
from airfoilShape import airfoilNACA16
from airfoilShape import airfoilBENZING
from airfoilShape import airfoilBICONVEX
from airfoilShape import airfoilMY_FILE
from airfoilShape import getNpt
import copy
import math as mt


def multiGeometry(nairfoils, n, params, alpha, dist, crel, typ, TE, ref_airfoil):
    # Multigeometry generation
    if ref_airfoil == "DEFAULT":
        ref_airfoil = np.argmax(crel) + 1

    crel_true = copy.copy(crel)  # copy of relative respective chords
    crel_true_y = copy.copy(crel)
    crel_true_y[0] = 0

    n_points = n

    index_Above1, index_Below1, index_Above2, index_Below2 = [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [None] * nairfoils
    x, y = [None] * nairfoils, [None] * nairfoils
    pointer, index_sx = [None] * nairfoils, [None] * nairfoils
    TE_len = [None] * nairfoils
    # coordinates x,y: blanked vectors

    for i in range(0, nairfoils):

        n = int(getNpt(n_points, crel[i]))  # n. of points of each airfoil, scaled by relative chord

        # get shape

        if typ[i] == "IGP":
            xcurr, ycurr, TE[i] = airfoilIGP(params[i], n, TE[i])  # x,y parameterized by IGP
        elif typ[i] == "NACA4":
            xcurr, ycurr = airfoilNACA4(params[i], n, TE[i])  # x,y parameterized using NACA 4-digit equation
        elif typ[i] == "NACA5":
            xcurr, ycurr = airfoilNACA5(params[i], n, TE[i])  # x,y parameterized using NACA 5-digit equation
        elif typ[i] == "NACA4_MOD":
            xcurr, ycurr = airfoilNACA4_MOD(params[i], n, TE[i])  # x,y parameterized using NACA 4-digit modified equation
        elif typ[i] == "NACA16":
            xcurr, ycurr = airfoilNACA16(params[i], n, TE[i])
        elif typ[i] == "BENZING":
            xcurr, ycurr, TE[i] = airfoilBENZING(params[i], n, TE[i], i)
        #elif typ[i] == "BICONVEX":
        #    xcurr, ycurr = airfoilBICONVEX(params[i], n, TE[i])
        elif typ[i] == ("MY_FILE_%d" % (i+1)):
            xcurr, ycurr, TE[i] = airfoilMY_FILE(i+1, TE[i])
        else:
            print('Error: parametrization not available for the airfoil number ' + str(i) +'.\n')


        xorig = copy.copy(xcurr.tolist())

        index_sx[i] = np.argmin(xcurr)
        index_Below1[i] = np.argmin(abs(np.asarray([xcurr[k] for k in range(0, index_sx[i] + 1)]) - 0.1))
        index_Above1[i] = index_sx[i] + np.argmin(
            abs(np.asarray([xcurr[k] for k in range(index_sx[i], len(xcurr))]) - xorig[index_Below1[i]]))
        index_Below2[i] = np.argmin(abs(np.asarray([xcurr[k] for k in range(0, index_sx[i] + 1)]) - 0.8))
        index_Above2[i] = index_sx[i] + np.argmin(
            abs(np.asarray([xcurr[k] for k in range(index_sx[i], len(xcurr))]) - xorig[index_Below2[i]]))

        if i > 0:
            xcurr = crel[i] * xcurr  # scaling x,y by respective relative chord
            ycurr = crel[i] * ycurr

            alpha_rel = alpha[i] - alpha[ref_airfoil - 1]  # relative AOA

            crel_true[i] = crel[i] * np.cos(alpha_rel * np.pi / 180)  # rotation of chord: x component
            crel_true_y[i] = crel_true_y[i] * np.sin(-alpha_rel * np.pi / 180)  # rotation of chord: y component

            # R: rotation matrix
            R = np.array([[np.cos(alpha_rel * np.pi / 180), np.sin(alpha_rel * np.pi / 180)], [-np.sin(alpha_rel * np.pi / 180), np.cos(alpha_rel * np.pi / 180)]])

            coord_mat = np.array([xcurr, ycurr])
            coord_mat = R.dot(coord_mat)  # new x,y coordinates: rotated

            xcurr = (coord_mat[0] + sum(crel_true[:i]) + sum([dist[j][0] for j in range(0, i + 1)])).tolist()
            ycurr = (coord_mat[1] + sum(crel_true_y[:i]) + sum([dist[j][1] for j in range(0, i + 1)])).tolist()

        x[i] = xcurr
        y[i] = ycurr

        if TE[i] == 'Y':
            TE_len[i] = mt.sqrt((x[i][len(x[i]) - 1] - x[i][0]) ** 2 + (y[i][len(x[i]) - 1] - y[i][0]) ** 2)
        else:
            pass

    return x, y, index_Above1, index_Below1, index_Above2, index_Below2, TE_len, ref_airfoil
