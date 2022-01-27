import numpy as np
import math as mt
from pts_gen import airfoil_pts
from pts_gen import grid_pts
from pts_gen import ext_pts_VISCOUS
from pts_gen import ext_pts_EULER
import copy


def PointsInCircum(r, n=100):
    x = [mt.cos(2 * mt.pi / n * x) * r for x in range(0, n + 1)]
    y = [mt.sin(2 * mt.pi / n * x) * r for x in range(0, n + 1)]
    return x, y


def mesh_gen_VISCOUS(x, y, nairfoils, index_Above1, index_Below1, index_Above2, index_Below2, alpha, dist, crel, TE, BL_thickness, norm_nodes, s, progr, airfoil_nodes, ellipse_nodes, h, nodes, norm_nodes_TE, TE_len, ellipse_dimension, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, ref_airfoil, farfield_size):
    # 1) BASIC VARIABLES DEFINITION AND INTRODUCTION
    # 2) MAIN POINTS OF EACH AIRFOIL AND STRUCTURED REGION
    #   2.1) AIRFOILS' POINTS
    #   2.2) STRUCTURED REGIONS' POINTS
    # 3) ELLIPSE
    #   3.1) DEFINITION AND MAIN POINTS
    #   3.2) ELLIPSE NODES CALCULUS
    # 4) FILE .GEO WRITING
    #   4.1) BASIC FARFIELD POINTS AND REFERENCES
    #   4.3) ELLIPSE'S POINTS, LINES AND SURFACE
    #   4.4) AIRFOILS' SPLINES, SURFACES AND STRUCTURED REGIONS
    #   4.5) FARFIELD
    #   4.6) ELLIPSE SURFACE DEFINITION & CONCLUSION OF REFINING POINTS and AIRFOILS' SURFACES
    #   4.7) REFINING POINTS (WAKE)
    #   4.8) PHYSICAL SURFACE AND LINES

    # ========================================================= 1) BASIC VARIABLES DEFINITION AND INTRODUCTION
    R = farfield_size  # radius farfield
    H = 1  # element size - farfield
    N_wake = 50
    thetaLE = mt.pi / 8
    thetaTE = mt.pi / 12

    # Empty variables definition
    old_gamma = [None] * (nairfoils + 1)
    N = [None] * nairfoils
    a1, x_wake1 = [None] * N_wake, [None] * N_wake
    a2, x_wake2 = [None] * N_wake, [None] * N_wake
    hh = [None] * N_wake
    x1, y1, x2, y2 = [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [None] * nairfoils
    m1, m2 = [None] * nairfoils, [None] * nairfoils
    bump = [None] * nairfoils
    indicator = [0] * nairfoils
    x_orig25, y_orig25, x_orig75, y_orig75 = [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [
        None] * nairfoils
    x_TE1, y_TE1, x_TE2, y_TE2, x_TE3, y_TE3, x_TE4, y_TE4, x_TE5, y_TE5, x_TE6, y_TE6, x_TE7, y_TE7, x_TE8, y_TE8 = [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils
    d = [None] * (nairfoils)

    # ========================================================= 2) MAIN POINTS OF EACH AIRFOIL AND STRUCTURED REGION
    # 2.1) AIRFOILS' POINTS
    y_ymin, y_ymax, x_sx, x_dx, index_ymin, index_ymax, index_sx, index_dx, curveLE, curveU, curveL, curveTE_U, curveTE_L, PSX, PDX, lengthLE, lengthUPPER, lengthLOWER, lengthTE_U, lengthTE_L = airfoil_pts(
        x, y, nairfoils, index_Above1, index_Below1, index_Above2, index_Below2, crel, airfoil_nodes)

    # 2.2) STRUCTURED REGIONS' POINTS
    xc, yc, theta, lengthLE_C, lengthUPPER_C, lengthLOWER_C, lengthTE_U_C, lengthTE_L_C, pointer = grid_pts(x, y,
                                                                                                            nairfoils,
                                                                                                            BL_thickness,
                                                                                                            index_sx,
                                                                                                            index_Above1,
                                                                                                            index_Below1,
                                                                                                            index_Above2,
                                                                                                            index_Below2,
                                                                                                            TE)

    for i in range(0, nairfoils):
        bump[i] = 5 + 218.72643 * ((lengthLE_C[i] - lengthLE[i]) ** 2) - 49.99787 * (lengthLE_C[i] - lengthLE[i])

    # ========================================================= 3) ELLIPSE
    # 3.1) DEFINITION AND MAIN POINTS
    centre_x, centre_y, x_ellipse, y_ellipse, sel_data, tau, ecc = ext_pts_VISCOUS(x, y, xc, yc, nairfoils, index_sx, PSX, PDX,
                                                                           crel, ellipse_dimension)

    mindiff = 1000000
    ellipse_area = mt.pi * sel_data.a * sel_data.b
    print('Ellipse Area: %.3f' % ellipse_area)

    # 0.25 chord of each airfoil
    for i in range(0, nairfoils):
        x_orig25[i] = (x[ref_airfoil - 1][i] - x[ref_airfoil - 1][index_sx[i]]) * 0.25 + x[ref_airfoil - 1][index_sx[i]]
        y_orig25[i] = (y[ref_airfoil - 1][i] - y[ref_airfoil - 1][index_sx[i]]) * 0.25 + y[ref_airfoil - 1][index_sx[i]]
        x_orig75[i] = (x[ref_airfoil - 1][i] - x[ref_airfoil - 1][index_sx[i]]) * 0.75 + x[ref_airfoil - 1][index_sx[i]]
        y_orig75[i] = (y[ref_airfoil - 1][i] - y[ref_airfoil - 1][index_sx[i]]) * 0.75 + y[ref_airfoil - 1][index_sx[i]]

    for i in range(0, len(x_ellipse)):
        diff = abs(y_ellipse[i] - y_orig25[ref_airfoil - 1] + mt.tan(alpha[ref_airfoil - 1] * mt.pi / 180) * (
                    x_orig25[ref_airfoil - 1] - x_ellipse[i]))

        if diff <= mindiff and i <= np.argmin(y_ellipse) and i >= np.argmax(y_ellipse) and x_ellipse[i] <= x_orig25[
            ref_airfoil - 1]:
            mindiff = diff
            old_gamma[0] = [i, 'LE']
    # ======================================
    indw1 = 0
    indw2 = 0
    for el in range(0, nairfoils):
        mindiff2 = 1000000

        if nairfoils == 1 or el == ref_airfoil - 1:
            q = y_orig25[el] - mt.tan(alpha[el] * mt.pi / 180) * x_orig25[el]

            for i in range(0, len(x_ellipse)):

                if alpha[el] >= 0:
                    if x_ellipse[i] > (x[el][index_sx[el]] + 0.75 * abs((x[el][0] - x[el][index_sx[el]]))) and \
                            y_ellipse[i] >= max(yc[el][0], yc[el][len(yc[el]) - 1]):
                        diff2 = abs(y_ellipse[i] - (mt.tan(alpha[el] * mt.pi / 180) * x_ellipse[i] + q))
                        if diff2 <= mindiff2:
                            mindiff2 = diff2
                            old_gamma[el + 1] = [i, 'TE']
                            indw1 = i
                else:
                    if x_ellipse[i] > (x[el][index_sx[el]] + 0.75 * abs((x[el][0] - x[el][index_sx[el]]))) and \
                            y_ellipse[i] <= min(yc[el][0], yc[el][len(yc[el]) - 1]):
                        diff2 = abs(y_ellipse[i] - (mt.tan(alpha[el] * mt.pi / 180) * x_ellipse[i] + q))
                        if diff2 <= mindiff2:
                            mindiff2 = diff2
                            old_gamma[el + 1] = [i, 'TE']
                            indw1 = i

            if nairfoils == 1:
                mindiff3 = 1000000
                old_gamma_new = 0
                q1 = y_orig75[el] - mt.tan(- alpha[el] * mt.pi / 180) * x_orig75[el]

                for i in range(0, len(x_ellipse)):

                    if alpha[el] >= 0:
                        if x_ellipse[i] >= x[el][0] and y_ellipse[i] <= min(yc[el][0], yc[el][len(yc[el]) - 1]):
                            diff3 = abs(y_ellipse[i] - (mt.tan(- alpha[el] * mt.pi / 180) * x_ellipse[i] + q1))

                            if diff3 <= mindiff3:
                                mindiff3 = diff3
                                old_gamma_new = [i, 'TE']
                                indw2 = i
                    else:
                        if x_ellipse[i] >= x[el][0] and y_ellipse[i] >= max(yc[el][0], yc[el][len(yc[el]) - 1]):
                            diff3 = abs(y_ellipse[i] - (mt.tan(- alpha[el] * mt.pi / 180) * x_ellipse[i] + q1))

                            if diff3 <= mindiff3:
                                mindiff3 = diff3
                                old_gamma_new = [i, 'TE']
                                indw2 = i

                old_gamma.append(old_gamma_new)


        elif el != PDX and el != ref_airfoil - 1:
            mindiff2 = 1000000
            q = y_orig75[el] - mt.tan((- alpha[el] + alpha[ref_airfoil - 1]) * mt.pi / 180) * x_orig75[el]

            for i in range(0, len(x_ellipse)):

                if x_ellipse[i] > x_orig75[el] and (i <= np.argmax(y_ellipse) or i >= np.argmin(y_ellipse)):

                    diff2 = abs(y_ellipse[i] - (
                                mt.tan((- alpha[el] + alpha[ref_airfoil - 1]) * mt.pi / 180) * x_ellipse[i] + q))

                    if diff2 <= mindiff2:
                        mindiff2 = diff2
                        old_gamma[el + 1] = [i, 'TE']

        else:
            mindiff2 = 1000000

            q1 = y_orig75[el] - mt.tan((- alpha[el] + alpha[ref_airfoil - 1]) * 1.5 * mt.pi / 180) * x_orig75[el]

            for i in range(0, len(x_ellipse)):

                if (alpha[el] - alpha[ref_airfoil - 1]) > 0:
                    if x_ellipse[i] >= x[el][0] and y_ellipse[i] <= min(yc[el]):

                        diff2 = abs(y_ellipse[i] - (
                                    mt.tan((-alpha[el] + alpha[ref_airfoil - 1]) * 1.5 * mt.pi / 180) * x_ellipse[
                                i] + q1))

                        if diff2 <= mindiff2:
                            mindiff2 = diff2
                            old_gamma[el + 1] = [i, 'TE']
                            indw2 = i

                elif (alpha[el] - alpha[ref_airfoil - 1]) == 0:
                    if 0 <= indw1 <= 180:
                        if x_ellipse[i] >= x[el][0] and y_ellipse[i] <= min(yc[el]):

                            diff2 = abs(y_ellipse[i] - (
                                        mt.tan((-alpha[el] + alpha[ref_airfoil - 1]) * 1.5 * mt.pi / 180) * x_ellipse[
                                    i] + q1))

                            if diff2 <= mindiff2:
                                mindiff2 = diff2
                                old_gamma[el + 1] = [i, 'TE']
                                indw2 = i
                    else:
                        if x_ellipse[i] >= x[el][0] and y_ellipse[i] >= max(yc[el]):
                            diff2 = abs(y_ellipse[i] - (
                                    mt.tan((-alpha[el] + alpha[ref_airfoil - 1]) * 1.5 * mt.pi / 180) * x_ellipse[
                                i] + q1))

                            if diff2 <= mindiff2:
                                mindiff2 = diff2
                                old_gamma[el + 1] = [i, 'TE']
                                indw2 = i

                else:
                    if x_ellipse[i] >= x[el][0] and y_ellipse[i] >= max(yc[el]):

                        diff2 = abs(y_ellipse[i] - (
                                mt.tan((-alpha[el] + alpha[ref_airfoil - 1]) * 1.5 * mt.pi / 180) * x_ellipse[
                            i] + q1))

                        if diff2 <= mindiff2:
                            mindiff2 = diff2
                            old_gamma[el + 1] = [i, 'TE']
                            indw2 = i

    gamma = sorted(old_gamma, key=lambda x: x[0])
    dummy_gamma = copy.copy(gamma)
    # Remove duplicates
    for g in range(0, len(gamma) - 1):
        if dummy_gamma[g][0] == dummy_gamma[g + 1][0]:
            gamma.remove(gamma[g + 1])
        else:
            pass

    # 3.2) ELLIPSE NODES CALCULUS
    ellnodes = [None] * (len(gamma) - 1)
    for i in range(0, len(gamma) - 1):
        if gamma[i + 1][1] == 'LE' or gamma[i][1] == 'LE':

            ellnodes[i] = int(abs(gamma[i + 1][0] - gamma[i][0]) * ellipse_nodes[0])

        else:

            ellnodes[i] = int(abs(gamma[i + 1][0] - gamma[i][0]) * ellipse_nodes[1])

    if len(gamma) == 2:
        ellnodes.append(int(abs(360 - gamma[len(gamma) - 1][0] + gamma[0][0]) * ellipse_nodes[0]))

    else:
        ellnodes.append(int(abs(360 - gamma[len(gamma) - 1][0] + gamma[0][0]) * ellipse_nodes[1]))

    # ========================================================= 4) FILE .GEO WRITING
    # 4.1) BASIC FARFIELD POINTS AND REFERENCES
    fid = open('Test.geo', 'w+')
    fid.write('h = %.10f; \n' % h)
    fid.write('H = %.10f; \n' % H)
    fid.write('R = %.3f; \n' % R)
    fid.write('pi = %.5f; \n' % mt.pi)
    fid.write('thetaLE = %.3f; \n' % thetaLE)
    fid.write('thetaTE = %.3f; \n\n' % thetaTE)
    fid.write('yLE = %.3f; \n' % y[PSX][index_ymax[PSX]])
    fid.write('xLE = %.3f; \n' % x[PSX][index_ymax[PSX]])
    fid.write('DxLE = (%.3f-yLE) * Tan(thetaLE);\n' % R)
    fid.write('yTE = %.3f; //TE point 1\n' % y[PDX][index_dx[PDX]])
    fid.write('xTE = %.3f; //TE point 1\n' % x[PDX][index_dx[PDX]])
    fid.write('DxTE = (%.3f-yTE) * Tan(thetaTE);\n\n' % R)
    fid.write('// farfield\n')
    fid.write('Point(99) = {0, 0, 0, h};\n')
    fid.write('Point(1) = {xLE-DxLE, %.3f, 0, H};\n' % R)
    fid.write('Point(2) = {xLE-DxLE, -%.3f, 0, H};\n' % R)
    fid.write('Point(3) = {%.3f, -%.3f, 0, H};\n' % (R, R))
    fid.write('Point(4) = {%.3f, %.3f, 0, H};\n\n' % (R, R))

    # 4.2) AIRFOILS' POINTS
    SUM = 0
    N[0] = copy.copy(
        SUM)  # N saves the number of points of the previous airfoils so that a continuous definition is possible.

    if nairfoils > 1:
        for p in range(1, nairfoils):
            SUM = SUM + len(x[p - 1]) + pointer[p - 1]
            N[p] = copy.copy(SUM)

    for k in range(0, nairfoils):

        for k1 in range(0, len(x[k])):
            fid.write('Point(100 + %d) = {%.12f,%.12f,0,h}; \n' % (N[k] + k1, x[k][k1], y[k][k1]))

        for k2 in range(0, len(xc[k])):
            fid.write('Point(10100 + %d) = {%.12f,%.12f,0,h}; \n' % (N[k] + k2, xc[k][k2], yc[k][k2]))

    # 4.3) ELLIPSE'S POINTS, LINES AND SURFACE
    for m in range(0, len(x_ellipse)):
        fid.write('Point(100000 + %d) = {%.12f,%.12f,0,h}; \n' % (m, x_ellipse[m], y_ellipse[m]))


    if len(gamma) == 2:
        fid.write(
            'Spline(5) = {(100000 + %d): (100000 + %d)};\n' % (gamma[0][0], gamma[1][0]))
        fid.write(
            'Spline(5 + 1) = {(100000 + %d):100359, 100000:(100000 + %d)};\n' % (gamma[1][0], gamma[0][0]))
        if ecc >= 0.01:
            fid.write('Transfinite Line{5} = %d Using Bump %.5f;\n' % (ellnodes[0], ecc / 4))
            fid.write('Transfinite Line{5 + 1} = %d Using Bump %.5f;\n' % (ellnodes[1], ecc / 4))
        else:
            fid.write('Transfinite Line{5} = %d Using Bump %.5f;\n' % (ellnodes[0], 0.0025))
            fid.write('Transfinite Line{5 + 1} = %d Using Bump %.5f;\n' % (ellnodes[1], 0.0025))
    else:
        for ell in range(0, len(gamma)):

            if ell < (len(gamma) - 1):

                if gamma[ell][1] == 'LE' or gamma[ell + 1][1] == 'LE':
                    fid.write(
                        'Spline(5 + %d) = {(100000 + %d): (100000 + %d)};\n' % (ell, gamma[ell][0], gamma[ell + 1][0]))

                    if ecc >= 0.01:
                        fid.write('Transfinite Line{5 + %d} = %d Using Bump %.5f;\n' % (ell, ellnodes[ell], ecc / 5))

                    else:
                        fid.write('Transfinite Line{5 + %d} = %d Using Bump %.5f;\n' % (ell, ellnodes[ell], 0.002))
                else:
                    fid.write(
                        'Spline(5 + %d) = {(100000 + %d): (100000 + %d)};\n' % (ell, gamma[ell][0], gamma[ell + 1][0]))
                    fid.write('Transfinite Line{5 + %d} = %d;\n' % (ell, ellnodes[ell]))

            else:

                if gamma[ell][1] == 'LE' or gamma[0][1] == 'LE':
                    fid.write('Spline(5 + %d) = {(100000 + %d):100359, 100000:(100000 + %d)};\n' % (
                        ell, gamma[ell][0], gamma[0][0]))

                    if ecc >= 0.01:
                        fid.write('Transfinite Line{5 + %d} = %d Using Bump %.5f;\n' % (ell, ellnodes[ell], ecc / 5))

                    else:
                        fid.write('Transfinite Line{5 + %d} = %d Using Bump %.5f;\n' % (ell, ellnodes[ell], 0.002))
                else:
                    fid.write('Spline(5 + %d) = {(100000 + %d):100359, 100000:(100000 + %d)};\n' % (
                        ell, gamma[ell][0], gamma[0][0]))
                    fid.write('Transfinite Line{5 + %d} = %d;\n' % (ell, ellnodes[ell]))

    fid.write('Curve Loop(2000) = {5, ')
    for j in range(1, len(gamma)):
        if j <= (len(gamma) - 2):
            fid.write('5 + %d, ' % j)
        else:
            fid.write('5 + %d};\n\n' % j)

    # 4.4) AIRFOILS' SPLINES, SURFACES AND STRUCTURED REGIONS
    for j in range(1, nairfoils + 1):
        i = j - 1

        fid.write('Spline(100*%d) = {100+%d:100+%d};\n' % (j, N[i], (N[i] + index_Below2[i])))
        fid.write(
            'Spline(100*%d + 1) = {100+%d:100+%d};\n' % (j, (N[i] + index_Below2[i]), (N[i] + index_Below1[i])))

        fid.write(
            'Spline(100*%d + 2) = {100+%d:100+%d};\n' % (j, (N[i] + index_Below1[i]), (N[i] + index_Above1[i])))

        fid.write(
            'Spline(100*%d + 3) = {100+%d:100+%d};\n' % (j, (N[i] + index_Above1[i]), (N[i] + index_Above2[i])))

        if TE[i] == 'Y':
            fid.write('Spline(100*%d + 4) = {100+%d:100+%d};\n' % (j, N[i] + index_Above2[i], N[i] + len(x[i]) - 1))
            fid.write('Line(100*%d + 5) = {100+%d, 100+%d};\n' % (j, N[i] + len(x[i]) - 1, N[i]))
            fid.write('Line(1000*%d + 11) = {100+%d,10100+%d};\n' % (j, N[i] + len(xc[i]) - 1, N[i] + len(xc[i]) - 1))

        else:
            fid.write('Spline(100*%d + 4) = {100+%d:100+%d, 100+%d};\n' % (
                j, N[i] + index_Above2[i], N[i] + len(x[i]) - 1, N[i]))
            fid.write('Line(1000*%d + 11) = {100+%d,10100+%d};\n' % (j, N[i], N[i] + len(xc[i]) - 1))

        fid.write('Spline(1000*%d) = {10100+%d:10100+%d};\n' % (j, N[i], (N[i] + index_Below2[i])))
        fid.write(
            'Spline(1000*%d + 1) = {10100+%d:10100+%d};\n' % (j, (N[i] + index_Below2[i]), (N[i] + index_Below1[i])))
        fid.write(
            'Spline(1000*%d + 2) = {10100+%d:10100+%d};\n' % (j, (N[i] + index_Below1[i]), (N[i] + index_Above1[i])))
        fid.write(
            'Spline(1000*%d + 3) = {10100+%d:10100+%d};\n' % (j, (N[i] + index_Above1[i]), (N[i] + index_Above2[i])))

        fid.write('Spline(1000*%d + 4) = {10100+%d:10100+%d};\n' % (j, N[i] + index_Above2[i], N[i] + len(xc[i]) - 1))

        fid.write('Line(1000*%d + 6) = {100+%d,10100+%d};\n' % (j, N[i], N[i]))
        fid.write(
            'Line(1000*%d + 7) = {100+%d,10100+%d};\n' % (j, (N[i] + index_Below2[i]), (N[i] + index_Below2[i])))

        fid.write(
            'Line(1000*%d + 8) = {100+%d,10100+%d};\n' % (j, (N[i] + index_Below1[i]), (N[i] + index_Below1[i])))

        fid.write(
            'Line(1000*%d + 9) = {100+%d,10100+%d};\n' % (j, (N[i] + index_Above1[i]), (N[i] + index_Above1[i])))

        fid.write(
            'Line(1000*%d + 10) = {100+%d,10100+%d};\n' % (j, (N[i] + index_Above2[i]), (N[i] + index_Above2[i])))

        fid.write('Transfinite Line{1000*%d, 100*%d} = %d Using Progression 1.015;\n' % (j, j, curveTE_L[i]))
        fid.write('Transfinite Line{1000*%d + 1, 100*%d + 1} = %d;\n' % (j, j, curveL[i]))

        fid.write('Transfinite Line{100*%d + 2} = %d Using Bump 5;\n' % (j, curveLE[i]))
        fid.write('Transfinite Line{1000*%d + 2} = %d Using Bump %.12f;\n' % (j, curveLE[i], bump[i]))

        fid.write('Transfinite Line{1000*%d + 3, 100*%d + 3} = %d;\n' % (j, j, curveU[i]))
        fid.write('Transfinite Line{+1000*%d + 4, 100*%d + 4} = %d Using Progression 0.985;\n' % (j, j, curveTE_U[i]))

        fid.write(
            'Transfinite Line{1000*%d + 6, 1000*%d + 7, 1000*%d + 8, 1000*%d + 9, 1000*%d + 10, 1000*%d + 11} = %d Using Progression %f;\n' % (
                j, j, j, j, j, j, norm_nodes[i], progr[i]))

        if TE[i] == 'Y':
            fid.write('Transfinite Line{100*%d + 5} = %d;\n' % (j, int(TE_len[i] * 20000 * nodes / crel[i])))
        else:
            pass

        fid.write('Curve Loop(10*%d) = {1000*%d, -1000*%d - 7, -100*%d, 1000*%d + 6};\n' % (j, j, j, j, j))
        fid.write('Curve Loop(10*%d + 1) = {100*%d + 1, 1000*%d + 8, -1000*%d - 1, -1000*%d - 7};\n' % (j, j, j, j, j))
        fid.write('Curve Loop(10*%d + 2) = {100*%d + 2, 1000*%d + 9, -1000*%d - 2, -1000*%d - 8};\n' % (j, j, j, j, j))
        fid.write('Curve Loop(10*%d + 3) = {100*%d + 3, 1000*%d + 10, -1000*%d - 3, -1000*%d - 9};\n' % (j, j, j, j, j))
        fid.write(
            'Curve Loop(10*%d + 4) = {100*%d + 4, 1000*%d + 11, -1000*%d - 4, -1000*%d - 10};\n' % (j, j, j, j, j))

        fid.write('Plane Surface(10*%d) = {10*%d};\n' % (j, j))
        fid.write('Plane Surface(10*%d + 1) = {10*%d + 1};\n' % (j, j))
        fid.write('Plane Surface(10*%d + 2) = {10*%d + 2};\n' % (j, j))
        fid.write('Plane Surface(10*%d + 3) = {10*%d + 3};\n' % (j, j))
        fid.write('Plane Surface(10*%d + 4) = {10*%d + 4};\n' % (j, j))

        fid.write('Transfinite Surface{10*%d};\n' % j)
        fid.write('Recombine Surface{10*%d};\n' % j)
        fid.write('Transfinite Surface{10*%d + 1};\n' % j)
        fid.write('Recombine Surface{10*%d + 1};\n' % j)
        fid.write('Transfinite Surface{10*%d + 2};\n' % j)
        fid.write('Recombine Surface{10*%d + 2};\n' % j)
        fid.write('Transfinite Surface{10*%d + 3};\n' % j)
        fid.write('Recombine Surface{10*%d + 3};\n' % j)
        fid.write('Transfinite Surface{10*%d + 4};\n' % j)
        fid.write('Recombine Surface{10*%d + 4};\n' % j)

    # 4.5) FARFIELD
    fid.write('Circle(1) = {1,99,2};\n')
    fid.write('Line(2) = {2,3};\n')
    fid.write('Line(3) = {3,4};\n')
    fid.write('Line(4) = {4,1};\n')

    nodes_farfield = [int(30 / max(crel) * nodes), int(25 / max(crel) * nodes), int(50 / max(crel) * nodes)]
    fid.write('Transfinite Line{1} = %d Using Bump 2;\n' % nodes_farfield[0])
    fid.write('Transfinite Line{2,-4} = %d Using Progression 0.995;\n' % nodes_farfield[1])
    fid.write('Transfinite Line{3} = %d Using Bump 3.5;\n\n' % nodes_farfield[2])

    fid.write('Curve Loop(1000) = {1,2,3,4};\n')

    fid.write('Plane Surface(1000) = {1000, 2000}; \n\n')

    # 4.6) WAKE
    if nairfoils == 1:
        angle_wake1 = alpha[ref_airfoil - 1]
        angle_wake2 = alpha[ref_airfoil - 1] * 0.5
    else:
        angle_wake1 = alpha[ref_airfoil - 1]
        angle_wake2 = 0
        sum_crel = 0

        for i in range(0, len(alpha)):
            angle_wake2 = angle_wake2 + (-alpha[i]) * crel[i]
            sum_crel = sum_crel + crel[i]
        angle_wake2 = angle_wake2 / sum_crel

    q_wake1 = y_ellipse[indw1] - mt.tan(angle_wake1 * mt.pi / 180) * x_ellipse[indw1]
    q_wake2 = y_ellipse[indw2] - mt.tan(angle_wake2 * mt.pi / 180) * x_ellipse[indw2]

    if sel_data.a >= sel_data.b:
        wake_len = wake_length * sel_data.a
    else:
        wake_len = wake_length * sel_data.b
    a1[0] = wake_len * (1.2 - 1) / ((1.2 ** N_wake) - 1)
    a2[0] = wake_len * (1.2 - 1) / ((1.2 ** N_wake) - 1)
    x_wake1[0] = x_ellipse[indw1] + a1[0]
    x_wake2[0] = x_ellipse[indw2] + a2[0]
    hh[0] = h

    for w in range(1, N_wake):
        a1[w] = a1[0] * (1.21 ** w)
        x_wake1[w] = x_wake1[w - 1] + a1[w]
        a2[w] = a2[0] * (1.2 ** w)
        x_wake2[w] = x_wake2[w - 1] + a2[w]
        if w <= 20:
            hh[w] = hh[0]
        elif (hh[0] * (wake_progr ** (w - 20))) < 0.5:
            hh[w] = hh[0] * (wake_progr ** (w - 20))
        else:
            hh[w] = 0.5

    # 4.6) ELLIPSE SURFACE DEFINITION & CONCLUSION OF REFINING POINTS and AIRFOILS' SURFACES
    # Its definition was possible only after the definition of each airfoil curve loop.

    for j in range(1, nairfoils + 1):
        if TE[j - 1] == 'Y':
            if abs(x[j - 1][0] - x[j - 1][len(x[j - 1]) - 1]) <= 1e-4:
                m1[j - 1] = 0
            else:
                m1[j - 1] = - 1 / ((y[j - 1][0] - y[j - 1][len(y[j - 1]) - 1]) / (
                        x[j - 1][0] - x[j - 1][len(x[j - 1]) - 1]))
            x1[j - 1] = x[j - 1][0] + BL_thickness[j - 1] * 0.01 * mt.cos(mt.atan(m1[j - 1]))
            y1[j - 1] = y[j - 1][0] + BL_thickness[j - 1] * 0.01 * mt.sin(mt.atan(m1[j - 1]))
            x2[j - 1] = x[j - 1][len(x[j - 1]) - 1] + BL_thickness[j - 1] * 0.01 * mt.cos(mt.atan(m1[j - 1]))
            y2[j - 1] = y[j - 1][len(x[j - 1]) - 1] + BL_thickness[j - 1] * 0.01 * mt.sin(mt.atan(m1[j - 1]))
            fid.write('Point(10*%d) = {%.12f,%.12f,0,h/4}; \n' % (j, x1[j - 1], y1[j - 1]))
            fid.write('Point(10*%d + 1) = {%.12f,%.12f,0,h/4}; \n' % (j, x2[j - 1], y2[j - 1]))

            fid.write('Line(100*%d + 6) = {100+%d, 10*%d};\n' % (j, N[j - 1], j))
            fid.write('Line(100*%d + 7) = {10*%d, 10*%d + 1};\n' % (j, j, j))
            fid.write('Line(100*%d + 8) = {10*%d + 1, 100 + %d};\n' % (j, j, N[j - 1] + len(x[j - 1]) - 1))
            fid.write('Transfinite Line{100*%d + 7} = %d;\n' % (
                j, int(TE_len[j - 1] * 20000 * nodes / crel[j - 1])))
            fid.write('Transfinite Line{100*%d + 6, -100*%d - 8} = %d Using Progression %f;\n' % (
                j, j, norm_nodes_TE[j - 1], progr[j-1]))
            fid.write('Curve Loop(10*%d + 5) = {100*%d + 5, 100*%d + 6, 100*%d + 7, 100*%d + 8};\n' % (j, j, j, j, j))
            fid.write('Plane Surface(10*%d + 5) = {10*%d + 5};\n' % (j, j))
            fid.write('Transfinite Surface{10*%d + 5};\n' % j)
            fid.write('Recombine Surface{10*%d + 5};\n' % j)

            fid.write('Curve Loop(100*%d) = {100*%d, 100*%d + 1, 100*%d + 2, 100*%d + 3, 100*%d + 4, 100*%d + 5};\n' % (
                j, j, j, j, j, j, j))
            fid.write(
                'Curve Loop(100*%d + 1) = {1000*%d + 6, 1000*%d, 1000*%d + 1, 1000*%d + 2, 1000*%d + 3, 1000*%d + 4, -1000*%d - 11, -100*%d - 8, -100*%d - 7, -100*%d - 6};\n' % (
                    j, j, j, j, j, j, j, j, j, j, j))

        else:
            m1[j - 1] = [None] * 2

            m1[j - 1][0] = (y[j - 1][0] - y[j - 1][1]) / (x[j - 1][0] - x[j - 1][1])
            m1[j - 1][1] = (y[j - 1][len(x[j - 1]) - 1] - y[j - 1][len(x[j - 1]) - 2]) / (
                    x[j - 1][len(x[j - 1]) - 1] - x[j - 1][len(x[j - 1]) - 2])
            x1[j - 1] = x[j - 1][0] + BL_thickness[j - 1] * 0.01 * mt.cos(mt.atan(m1[j - 1][0]))
            y1[j - 1] = y[j - 1][0] + BL_thickness[j - 1] * 0.01 * mt.sin(mt.atan(m1[j - 1][0]))
            x2[j - 1] = x[j - 1][0] + BL_thickness[j - 1] * 0.01 * mt.cos(mt.atan(m1[j - 1][1]))
            y2[j - 1] = y[j - 1][0] + BL_thickness[j - 1] * 0.01 * mt.sin(mt.atan(m1[j - 1][1]))
            fid.write('Point(10*%d) = {%.12f,%.12f,0,h/4}; \n' % (j, x1[j - 1], y1[j - 1]))
            fid.write('Point(10*%d + 1) = {%.12f,%.12f,0,h/4}; \n' % (j, x2[j - 1], y2[j - 1]))
            fid.write(
                'Curve Loop(100*%d) = {100*%d, 100*%d + 1, 100*%d + 2, 100*%d + 3, 100*%d + 4};\n' % (j, j, j, j, j, j))
            fid.write(
                'Curve Loop(100*%d + 1) = {1000*%d + 6, 1000*%d, 1000*%d + 1, 1000*%d + 2, 1000*%d + 3, 1000*%d + 4, -1000*%d - 11};\n' % (
                    j, j, j, j, j, j, j, j))

    fid.write('Plane Surface(2000) = {2000')
    for j in range(1, nairfoils + 1):
        if nairfoils == 1:
            fid.write(', 101}; \n\n')
        elif j != nairfoils:
            fid.write(', 100*%d + 1' % j)
        else:
            fid.write(', 100*%d + 1};\n\n' % j)

    for i in range(0, nairfoils):
        if TE[i] == 'Y':
            pass
        else:
            fid.write('Point{10*%d} In Surface {2000}; \n' % (i + 1))
            fid.write('Point{10*%d + 1} In Surface {2000}; \n' % (i + 1))

    for j in range(0, len(gamma)):
        fid.write('Point{100000 + %d} In Surface {1000}; \n' % gamma[j][0])
        fid.write('Point{100000 + %d} In Surface {2000}; \n' % gamma[j][0])

    # 4.7) REFINING POINTS (WAKE)
    for k in range(0, nairfoils):

        if nairfoils == 1 or k == (nairfoils - 1):
            mindiff4 = 1000000

            if nairfoils == 1:

                if TE[k] == 'Y':
                    q2 = (y[k][0] + y[k][len(x) - 1]) / 2 - mt.tan(alpha[k] * mt.pi / 180) * (
                            x[k][0] + x[k][len(x) - 1]) / 2
                else:
                    q2 = y[k][0] - mt.tan(alpha[k] * mt.pi / 180) * x[k][0]

                for i in range(0, len(x_ellipse)):
                    if x_ellipse[i] >= x[k][0]:
                        diff4 = abs(y_ellipse[i] - (mt.tan(alpha[k] * mt.pi / 180) * x_ellipse[i] + q2))

                        if diff4 <= mindiff4:
                            mindiff4 = diff4
                            indicator[k] = i

            else:
                if TE[k] == 'Y':
                    q2 = (y[k][0] + y[k][len(x) - 1]) / 2 - mt.tan(
                        (- alpha[k] + alpha[ref_airfoil - 1]) * mt.pi / 180) * (x[k][0] + x[k][len(x) - 1]) / 2
                else:
                    q2 = y[k][0] - mt.tan((- alpha[k] + alpha[ref_airfoil - 1]) * mt.pi / 180) * x[k][0]

                for i in range(0, len(x_ellipse)):
                    if x_ellipse[i] >= x[k][0]:
                        diff4 = abs(y_ellipse[i] - (
                                    mt.tan((-alpha[k] + alpha[ref_airfoil - 1]) * mt.pi / 180) * x_ellipse[i] + q2))

                        if diff4 <= mindiff4:
                            mindiff4 = diff4
                            indicator[k] = i

            if TE[k] == 'Y':
                if alpha[k] >= 0:
                    x_TE1[k] = (x[k][0] + x[k][1] + x[k][len(x[k]) - 1] + x[k][len(x[k]) - 2] * 2 + x_ellipse[
                        indicator[k]]) / 6
                    y_TE1[k] = (y[k][0] + y[k][1] + y[k][len(x[k]) - 1] + y[k][len(x[k]) - 2] * 2 + y_ellipse[
                        indicator[k]]) / 6
                else:
                    x_TE1[k] = (x[k][0] + x[k][1] * 2 + x[k][len(x[k]) - 1] + x[k][len(x[k]) - 2] + x_ellipse[
                        indicator[k]]) / 6
                    y_TE1[k] = (y[k][0] + y[k][1] * 2 + y[k][len(x[k]) - 1] + y[k][len(x[k]) - 2] + y_ellipse[
                        indicator[k]]) / 6
                x_TE2[k] = (x[k][0] + x[k][1] + x[k][len(x[k]) - 1] + x[k][len(x[k]) - 2] + x_ellipse[indicator[k]]) / 5
                y_TE2[k] = (y[k][0] + y[k][1] + y[k][len(x[k]) - 1] + y[k][len(x[k]) - 2] + y_ellipse[indicator[k]]) / 5
                x_TE3[k] = (x[k][0] + x[k][len(x[k]) - 1] + x_ellipse[indicator[k]]) / 3
                y_TE3[k] = (y[k][0] + y[k][len(x[k]) - 1] + y_ellipse[indicator[k]]) / 3
                x_TE4[k] = ((x[k][0] + x[k][len(x[k]) - 1]) / 2 + x_ellipse[indicator[k]]) / 2
                y_TE4[k] = ((y[k][0] + y[k][len(x[k]) - 1]) / 2 + y_ellipse[indicator[k]]) / 2
                x_TE5[k] = ((x[k][0] + x[k][len(x[k]) - 1]) / 2 + x_ellipse[indicator[k]] * 2) / 3
                y_TE5[k] = ((y[k][0] + y[k][len(x[k]) - 1]) / 2 + y_ellipse[indicator[k]] * 2) / 3
                d[k] = mt.sqrt(
                    ((x[k][0] + x[k][len(x[k]) - 1]) / 2 - x_ellipse[indicator[k]]) ** 2 + (
                            (y[k][0] + y[k][len(x[k]) - 1]) / 2 - y_ellipse[indicator[k]]) ** 2)
                x_TE6[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] / 5
                y_TE6[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2
                x_TE7[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] / 3
                y_TE7[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2
                x_TE8[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] * 2 / 3
                y_TE8[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2

            else:
                x_TE1[k] = (x[k][0] * 5 + x_ellipse[indicator[k]]) / 6
                y_TE1[k] = (y[k][0] * 5 + y_ellipse[indicator[k]]) / 6
                x_TE2[k] = (x[k][0] * 4 + x_ellipse[indicator[k]]) / 5
                y_TE2[k] = (y[k][0] * 4 + y_ellipse[indicator[k]]) / 5
                x_TE3[k] = (x[k][0] * 2 + x_ellipse[indicator[k]]) / 3
                y_TE3[k] = (y[k][0] * 2 + y_ellipse[indicator[k]]) / 3
                x_TE4[k] = (x[k][0] + x_ellipse[indicator[k]]) / 2
                y_TE4[k] = (y[k][0] + y_ellipse[indicator[k]]) / 2
                x_TE5[k] = (x[k][0] + x_ellipse[indicator[k]] * 2) / 3
                y_TE5[k] = (y[k][0] + y_ellipse[indicator[k]] * 2) / 3
                d[k] = mt.sqrt(
                    (x[k][0] - x_ellipse[indw2]) ** 2 + (
                            y[k][0] - y_ellipse[indw2]) ** 2)
                x_TE6[k] = x[k][0] + d[k] / 5
                y_TE6[k] = y[k][0]
                x_TE7[k] = x[k][0] + d[k] / 3
                y_TE7[k] = y[k][0]
                x_TE8[k] = x[k][0] + d[k] * 2 / 3
                y_TE8[k] = y[k][0]

        else:
            MINDIFF = 10000
            if xc[k + 1][index_sx[k + 1]] >= max(xc[k][0], xc[k][len(xc[k]) - 1]):
                X1 = x[k][np.argmax(x[k])]
                Y1 = y[k][np.argmax(x[k])]
                # X2 = xc[k + 1][index_sx[k + 1]]
                # Y2 = yc[k + 1][index_sx[k + 1]]

            elif xc[k + 1][index_sx[k + 1]] < max(xc[k][0], xc[k][len(xc[k]) - 1]) and yc[k + 1][index_sx[k + 1]] >= \
                    yc[k][len(xc[k]) - 1]:
                X1 = xc[k][len(xc[k]) - 1]
                Y1 = yc[k][len(xc[k]) - 1]
                # X2 = xc[k + 1][index_Below1[k + 1]]
                # Y2 = yc[k + 1][index_Below1[k + 1]]

            else:
                X1 = xc[k][0]
                Y1 = yc[k][0]
                # X2 = xc[k + 1][index_Above1[k + 1]]
                # Y2 = yc[k + 1][index_Above1[k + 1]]

            for i in range(0, len(xc[k + 1])):
                DIFF = mt.sqrt((X1 - xc[k + 1][i]) ** 2 + (Y1 - yc[k + 1][i]) ** 2)
                if DIFF <= MINDIFF:
                    MINDIFF = DIFF
                    ind_gap = i

            x_TE1[k] = (X1 * 5 + xc[k + 1][ind_gap]) / 6
            y_TE1[k] = (Y1 * 5 + yc[k + 1][ind_gap]) / 6

            if ind_gap >= 2:
                x_TE2[k] = (X1 * 4 + xc[k + 1][ind_gap - 1]) / 5
                y_TE2[k] = (Y1 * 4 + yc[k + 1][ind_gap - 1]) / 5
                x_TE3[k] = (X1 * 2 + xc[k + 1][ind_gap - 1]) / 3
                y_TE3[k] = (Y1 * 2 + yc[k + 1][ind_gap - 1]) / 3
                x_TE4[k] = (X1 + xc[k + 1][ind_gap - 2]) / 2
                y_TE4[k] = (Y1 + yc[k + 1][ind_gap - 2]) / 2
                x_TE5[k] = (X1 + xc[k + 1][ind_gap - 2] * 2) / 3
                y_TE5[k] = (Y1 + yc[k + 1][ind_gap - 2] * 2) / 3
            elif ind_gap >= 1:
                x_TE2[k] = (X1 * 4 + xc[k + 1][ind_gap - 1]) / 5
                y_TE2[k] = (Y1 * 4 + yc[k + 1][ind_gap - 1]) / 5
                x_TE3[k] = (X1 * 2 + xc[k + 1][ind_gap - 1]) / 3
                y_TE3[k] = (Y1 * 2 + yc[k + 1][ind_gap - 1]) / 3
                x_TE4[k] = (X1 + xc[k + 1][ind_gap - 1]) / 2
                y_TE4[k] = (Y1 + xc[k + 1][ind_gap - 1]) / 2
                x_TE5[k] = (X1 + yc[k + 1][ind_gap - 1] * 2) / 3
                y_TE5[k] = (Y1 + yc[k + 1][ind_gap - 1] * 2) / 3
            else:
                x_TE2[k] = (X1 * 4 + xc[k + 1][ind_gap]) / 5
                y_TE2[k] = (Y1 * 4 + yc[k + 1][ind_gap]) / 5
                x_TE3[k] = (X1 * 2 + xc[k + 1][ind_gap]) / 3
                y_TE3[k] = (Y1 * 2 + yc[k + 1][ind_gap]) / 3
                x_TE4[k] = (X1 + xc[k + 1][ind_gap]) / 2
                y_TE4[k] = (Y1 + yc[k + 1][ind_gap]) / 2
                x_TE5[k] = (X1 + xc[k + 1][ind_gap] * 2) / 3
                y_TE5[k] = (Y1 + yc[k + 1][ind_gap] * 2) / 3

            if TE[k] == 'Y':
                d[k] = mt.sqrt(
                    ((x[k][0] + x[k][len(x[k]) - 1]) / 2 - xc[k + 1][index_sx[k + 1]]) ** 2 + (
                            (y[k][0] + y[k][len(x[k]) - 1]) / 2 - yc[k + 1][index_sx[k + 1]]) ** 2)
                x_TE6[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] / 5
                y_TE6[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2
                x_TE7[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] / 3
                y_TE7[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2
                x_TE8[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] * 2 / 3
                y_TE8[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2

            else:
                d[k] = mt.sqrt(
                    (x[k][0] - xc[k + 1][index_sx[k + 1]]) ** 2 + (
                            y[k][0] - yc[k + 1][index_sx[k + 1]]) ** 2)
                x_TE6[k] = x[k][0] + d[k] / 5
                y_TE6[k] = y[k][0]
                x_TE7[k] = x[k][0] + d[k] / 3
                y_TE7[k] = y[k][0]
                x_TE8[k] = x[k][0] + d[k] * 2 / 3
                y_TE8[k] = y[k][0]

        fid.write('Point(10*%d + 2) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE1[k], y_TE1[k]))
        fid.write('Point{10*%d + 2} In Surface {2000}; \n' % (k + 1))
        fid.write('Point(10*%d + 3) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE2[k], y_TE2[k]))
        fid.write('Point{10*%d + 3} In Surface {2000}; \n' % (k + 1))
        fid.write('Point(10*%d + 4) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE3[k], y_TE3[k]))
        fid.write('Point{10*%d + 4} In Surface {2000}; \n' % (k + 1))
        fid.write('Point(10*%d + 5) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE4[k], y_TE4[k]))
        fid.write('Point{10*%d + 5} In Surface {2000}; \n' % (k + 1))
        fid.write('Point(10*%d + 6) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE5[k], y_TE5[k]))
        fid.write('Point{10*%d + 6} In Surface {2000}; \n' % (k + 1))
        fid.write('Point(10*%d + 7) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE6[k], y_TE6[k]))
        fid.write('Point{10*%d + 7} In Surface {2000}; \n' % (k + 1))
        fid.write('Point(10*%d + 8) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE7[k], y_TE7[k]))
        fid.write('Point{10*%d + 8} In Surface {2000}; \n' % (k + 1))
        fid.write('Point(10*%d + 9) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE8[k], y_TE8[k]))
        fid.write('Point{10*%d + 9} In Surface {2000}; \n' % (k + 1))

    ll = 0
    ind = 0
    while ind < N_wake:
        fid.write('Point(100000 + %d) = {%.12f,%.12f,0,%.12f}; \n' % (
            len(x_ellipse) + ll, x_wake1[ind], q_wake1 + mt.tan(angle_wake1 * mt.pi / 180) * x_wake1[ind], hh[ind]))
        fid.write('Point{100000 + %d} In Surface {1000}; \n' % (len(x_ellipse) + ll))
        ll = ll + 1
        fid.write('Point(100000 + %d) = {%.12f,%.12f,0,%.12f}; \n' % (
            len(x_ellipse) + ll, x_wake2[ind], q_wake2 + mt.tan(angle_wake2 * mt.pi / 180) * x_wake2[ind], hh[ind]))
        fid.write('Point{100000 + %d} In Surface {1000}; \n' % (len(x_ellipse) + ll))
        ll = ll + 1
        ind = ind + 1

    # 4.8) PHYSICAL SURFACE AND LINES
    fid.write('Physical Surface(1) = {1000, 2000')
    for j in range(1, nairfoils + 1):
        if nairfoils == 1:
            if TE[0] == 'Y':
                fid.write(', 10, 11, 12, 13, 14, 15};\n\n')
            else:
                fid.write(', 10, 11, 12, 13, 14};\n\n')
        elif j != nairfoils:
            if TE[j - 1] == 'Y':
                fid.write(', 10*%d, 10*%d + 1, 10*%d + 2, 10*%d + 3, 10*%d + 4, 10*%d + 5' % (j, j, j, j, j, j))
            else:
                fid.write(', 10*%d, 10*%d + 1, 10*%d + 2, 10*%d + 3, 10*%d + 4' % (j, j, j, j, j))
        else:
            if TE[j - 1] == 'Y':
                fid.write(', 10*%d, 10*%d + 1, 10*%d + 2, 10*%d + 3, 10*%d + 4, 10*%d + 5}; \n\n' % (j, j, j, j, j, j))
            else:
                fid.write(', 10*%d, 10*%d + 1, 10*%d + 2, 10*%d + 3, 10*%d + 4}; \n\n' % (j, j, j, j, j))

    for j in range(1, nairfoils + 1):
        if TE[j - 1] == 'N':
            fid.write('Physical Line("airfoil%d") = {100*%d, 100*%d + 1, 100*%d + 2, 100*%d + 3, 100*%d + 4};\n' % (
                j, j, j, j, j, j))
        else:
            fid.write(
                'Physical Line("airfoil%d") = {100*%d, 100*%d + 1, 100*%d + 2, 100*%d + 3, 100*%d + 4, 100*%d + 5};\n' % (
                    j, j, j, j, j, j, j))

    fid.write('Physical Line("farfield") = {1,2,3,4};\n')

    # File .su2 generation.
    fid.write('Mesh.Format = 42;\n')
    if Mesh_Algo == 6:
        fid.write('Mesh.Algorithm = 6;\n')
    elif Mesh_Algo == 'DEFAULT' or Mesh_Algo == 5:
        fid.write('Mesh.Algorithm = 5;\n')
    else:
        fid.write('Mesh.Algorithm = %d;\n' % Mesh_Algo)

    if external_pts == 'YES':
        x_circle_ext, y_circle_ext = PointsInCircum(min(max(sel_data.b * semicircle_dimension, sel_data.a * semicircle_dimension), 200 * max(crel)), 360)
        i = 90
        while i <= 270:
            fid.write('Point(100000 + %d) = {%.12f,%.12f,0,%f}; \n' % (len(x_ellipse) + ll, x_circle_ext[i], y_circle_ext[i], h * semicircle_elem_factor))
            fid.write('Point{100000 + %d} In Surface {1000}; \n' % (len(x_ellipse) + ll))
            ll = ll + 1
            i = i + 10

    else:
        pass

    fid.write('Mesh 2;\n')

    fid.close()

    return []


def mesh_gen_EULER(x, y, nairfoils, index_Above1, index_Below1, index_Above2, index_Below2, alpha, dist, crel, TE, airfoil_nodes, ellipse_nodes, h, nodes, TE_len, ellipse_dimension, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, ref_airfoil, farfield_size):
    # 1) BASIC VARIABLES DEFINITION AND INTRODUCTION
    # 2) MAIN POINTS OF EACH AIRFOIL AND STRUCTURED REGION
    #   2.1) AIRFOILS' POINTS
    #   2.2) STRUCTURED REGIONS' POINTS
    # 3) ELLIPSE
    #   3.1) DEFINITION AND MAIN POINTS
    #   3.2) ELLIPSE NODES CALCULUS
    # 4) FILE .GEO WRITING
    #   4.1) BASIC FARFIELD POINTS AND REFERENCES
    #   4.3) ELLIPSE'S POINTS, LINES AND SURFACE
    #   4.4) AIRFOILS' SPLINES, SURFACES AND STRUCTURED REGIONS
    #   4.5) FARFIELD
    #   4.6) ELLIPSE SURFACE DEFINITION & CONCLUSION OF REFINING POINTS and AIRFOILS' SURFACES
    #   4.7) REFINING POINTS (WAKE)
    #   4.8) PHYSICAL SURFACE AND LINES

    # ========================================================= 1) BASIC VARIABLES DEFINITION AND INTRODUCTION
    R = farfield_size  # radius farfield
    H = 1  # element size - farfield
    N_wake = 50
    thetaLE = mt.pi / 8
    thetaTE = mt.pi / 12

    # Empty variables definition
    old_gamma = [None] * (nairfoils + 1)
    N = [None] * nairfoils
    a1, x_wake1 = [None] * N_wake, [None] * N_wake
    a2, x_wake2 = [None] * N_wake, [None] * N_wake
    hh = [None] * N_wake
    x1, y1, x2, y2 = [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [None] * nairfoils
    m1, m2 = [None] * nairfoils, [None] * nairfoils
    bump = [None] * nairfoils
    indicator = [0] * nairfoils
    x_orig25, y_orig25, x_orig75, y_orig75 = [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [
        None] * nairfoils
    x_TE1, y_TE1, x_TE2, y_TE2, x_TE3, y_TE3, x_TE4, y_TE4, x_TE5, y_TE5, x_TE6, y_TE6, x_TE7, y_TE7, x_TE8, y_TE8 = [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils, [
                                                                                                                         None] * nairfoils
    d = [None] * (nairfoils)

    # ========================================================= 2) MAIN POINTS OF EACH AIRFOIL AND STRUCTURED REGION
    # 2.1) AIRFOILS' POINTS
    y_ymin, y_ymax, x_sx, x_dx, index_ymin, index_ymax, index_sx, index_dx, curveLE, curveU, curveL, curveTE_U, curveTE_L, PSX, PDX, lengthLE, lengthUPPER, lengthLOWER, lengthTE_U, lengthTE_L = airfoil_pts(
        x, y, nairfoils, index_Above1, index_Below1, index_Above2, index_Below2, crel, airfoil_nodes)

    # ========================================================= 3) ELLIPSE
    # 3.1) DEFINITION AND MAIN POINTS
    centre_x, centre_y, x_ellipse, y_ellipse, sel_data, tau, ecc = ext_pts_EULER(x, y, nairfoils, index_sx, PSX, PDX,
                                                                           crel, ellipse_dimension)

    mindiff = 1000000
    ellipse_area = mt.pi * sel_data.a * sel_data.b
    print('Ellipse Area: %.3f' % ellipse_area)

    # 0.25 chord of each airfoil
    for i in range(0, nairfoils):
        x_orig25[i] = (x[ref_airfoil - 1][i] - x[ref_airfoil - 1][index_sx[i]]) * 0.25 + x[ref_airfoil - 1][index_sx[i]]
        y_orig25[i] = (y[ref_airfoil - 1][i] - y[ref_airfoil - 1][index_sx[i]]) * 0.25 + y[ref_airfoil - 1][index_sx[i]]
        x_orig75[i] = (x[ref_airfoil - 1][i] - x[ref_airfoil - 1][index_sx[i]]) * 0.75 + x[ref_airfoil - 1][index_sx[i]]
        y_orig75[i] = (y[ref_airfoil - 1][i] - y[ref_airfoil - 1][index_sx[i]]) * 0.75 + y[ref_airfoil - 1][index_sx[i]]

    for i in range(0, len(x_ellipse)):
        diff = abs(y_ellipse[i] - y_orig25[ref_airfoil - 1] + mt.tan(alpha[ref_airfoil - 1] * mt.pi / 180) * (
                    x_orig25[ref_airfoil - 1] - x_ellipse[i]))

        if diff <= mindiff and i <= np.argmin(y_ellipse) and i >= np.argmax(y_ellipse) and x_ellipse[i] <= x_orig25[
            ref_airfoil - 1]:
            mindiff = diff
            old_gamma[0] = [i, 'LE']
    # ======================================
    indw1 = 0
    indw2 = 0
    for el in range(0, nairfoils):
        mindiff2 = 1000000

        if nairfoils == 1 or el == ref_airfoil - 1:
            q = y_orig25[el] - mt.tan(alpha[el] * mt.pi / 180) * x_orig25[el]

            for i in range(0, len(x_ellipse)):

                if alpha[el] >= 0:
                    if x_ellipse[i] > (x[el][index_sx[el]] + 0.75 * abs((x[el][0] - x[el][index_sx[el]]))):
                        diff2 = abs(y_ellipse[i] - (mt.tan(alpha[el] * mt.pi / 180) * x_ellipse[i] + q))
                        if diff2 <= mindiff2:
                            mindiff2 = diff2
                            old_gamma[el + 1] = [i, 'TE']
                            indw1 = i
                else:
                    if x_ellipse[i] > (x[el][index_sx[el]] + 0.75 * abs((x[el][0] - x[el][index_sx[el]]))):
                        diff2 = abs(y_ellipse[i] - (mt.tan(alpha[el] * mt.pi / 180) * x_ellipse[i] + q))
                        if diff2 <= mindiff2:
                            mindiff2 = diff2
                            old_gamma[el + 1] = [i, 'TE']
                            indw1 = i

            if nairfoils == 1:
                mindiff3 = 1000000
                old_gamma_new = 0
                q1 = y_orig75[el] - mt.tan(- alpha[el] * mt.pi / 180) * x_orig75[el]

                for i in range(0, len(x_ellipse)):

                    if alpha[el] >= 0:
                        if x_ellipse[i] >= x[el][0]:
                            diff3 = abs(y_ellipse[i] - (mt.tan(- alpha[el] * mt.pi / 180) * x_ellipse[i] + q1))

                            if diff3 <= mindiff3:
                                mindiff3 = diff3
                                old_gamma_new = [i, 'TE']
                                indw2 = i
                    else:
                        if x_ellipse[i] >= x[el][0]:
                            diff3 = abs(y_ellipse[i] - (mt.tan(- alpha[el] * mt.pi / 180) * x_ellipse[i] + q1))

                            if diff3 <= mindiff3:
                                mindiff3 = diff3
                                old_gamma_new = [i, 'TE']
                                indw2 = i

                old_gamma.append(old_gamma_new)


        elif el != PDX and el != ref_airfoil - 1:
            mindiff2 = 1000000
            q = y_orig75[el] - mt.tan((- alpha[el] + alpha[ref_airfoil - 1]) * mt.pi / 180) * x_orig75[el]

            for i in range(0, len(x_ellipse)):

                if x_ellipse[i] > x_orig75[el] and (i <= np.argmax(y_ellipse) or i >= np.argmin(y_ellipse)):

                    diff2 = abs(y_ellipse[i] - (
                                mt.tan((- alpha[el] + alpha[ref_airfoil - 1]) * mt.pi / 180) * x_ellipse[i] + q))

                    if diff2 <= mindiff2:
                        mindiff2 = diff2
                        old_gamma[el + 1] = [i, 'TE']

        else:
            mindiff2 = 1000000

            q1 = y_orig75[el] - mt.tan((- alpha[el] + alpha[ref_airfoil - 1]) * 1.5 * mt.pi / 180) * x_orig75[el]

            for i in range(0, len(x_ellipse)):

                if (alpha[el] - alpha[ref_airfoil - 1]) > 0:
                    if x_ellipse[i] >= x[el][0]:

                        diff2 = abs(y_ellipse[i] - (
                                    mt.tan((-alpha[el] + alpha[ref_airfoil - 1]) * 1.5 * mt.pi / 180) * x_ellipse[
                                i] + q1))

                        if diff2 <= mindiff2:
                            mindiff2 = diff2
                            old_gamma[el + 1] = [i, 'TE']
                            indw2 = i

                elif (alpha[el] - alpha[ref_airfoil - 1]) == 0:
                    if 0 <= indw1 <= 180:
                        if x_ellipse[i] >= x[el][0]:

                            diff2 = abs(y_ellipse[i] - (
                                        mt.tan((-alpha[el] + alpha[ref_airfoil - 1]) * 1.5 * mt.pi / 180) * x_ellipse[
                                    i] + q1))

                            if diff2 <= mindiff2:
                                mindiff2 = diff2
                                old_gamma[el + 1] = [i, 'TE']
                                indw2 = i
                    else:
                        if x_ellipse[i] >= x[el][0]:
                            diff2 = abs(y_ellipse[i] - (
                                    mt.tan((-alpha[el] + alpha[ref_airfoil - 1]) * 1.5 * mt.pi / 180) * x_ellipse[
                                i] + q1))

                            if diff2 <= mindiff2:
                                mindiff2 = diff2
                                old_gamma[el + 1] = [i, 'TE']
                                indw2 = i

                else:
                    if x_ellipse[i] >= x[el][0]:

                        diff2 = abs(y_ellipse[i] - (
                                mt.tan((-alpha[el] + alpha[ref_airfoil - 1]) * 1.5 * mt.pi / 180) * x_ellipse[
                            i] + q1))

                        if diff2 <= mindiff2:
                            mindiff2 = diff2
                            old_gamma[el + 1] = [i, 'TE']
                            indw2 = i

    gamma = sorted(old_gamma, key=lambda x: x[0])
    dummy_gamma = copy.copy(gamma)
    # Remove duplicates
    for g in range(0, len(gamma) - 1):
        if dummy_gamma[g][0] == dummy_gamma[g + 1][0]:
            gamma.remove(gamma[g + 1])
        else:
            pass

    # 3.2) ELLIPSE NODES CALCULUS
    ellnodes = [None] * (len(gamma) - 1)
    for i in range(0, len(gamma) - 1):
        if gamma[i + 1][1] == 'LE' or gamma[i][1] == 'LE':

            ellnodes[i] = int(abs(gamma[i + 1][0] - gamma[i][0]) * ellipse_nodes[0])

        else:

            ellnodes[i] = int(abs(gamma[i + 1][0] - gamma[i][0]) * ellipse_nodes[1])

    if len(gamma) == 2:
        ellnodes.append(int(abs(360 - gamma[len(gamma) - 1][0] + gamma[0][0]) * ellipse_nodes[0]))

    else:
        ellnodes.append(int(abs(360 - gamma[len(gamma) - 1][0] + gamma[0][0]) * ellipse_nodes[1]))

    # ========================================================= 4) FILE .GEO WRITING
    # 4.1) BASIC FARFIELD POINTS AND REFERENCES
    fid = open('Test.geo', 'w+')
    fid.write('h = %.10f; \n' % h)
    fid.write('H = %.10f; \n' % H)
    fid.write('R = %.3f; \n' % R)
    fid.write('pi = %.5f; \n' % mt.pi)
    fid.write('thetaLE = %.3f; \n' % thetaLE)
    fid.write('thetaTE = %.3f; \n\n' % thetaTE)
    fid.write('yLE = %.3f; \n' % y[PSX][index_ymax[PSX]])
    fid.write('xLE = %.3f; \n' % x[PSX][index_ymax[PSX]])
    fid.write('DxLE = (%.3f-yLE) * Tan(thetaLE);\n' % R)
    fid.write('yTE = %.3f; //TE point 1\n' % y[PDX][index_dx[PDX]])
    fid.write('xTE = %.3f; //TE point 1\n' % x[PDX][index_dx[PDX]])
    fid.write('DxTE = (%.3f-yTE) * Tan(thetaTE);\n\n' % R)
    fid.write('// farfield\n')
    fid.write('Point(99) = {0, 0, 0, h};\n')
    fid.write('Point(1) = {xLE-DxLE, %.3f, 0, H};\n' % R)
    fid.write('Point(2) = {xLE-DxLE, -%.3f, 0, H};\n' % R)
    fid.write('Point(3) = {%.3f, -%.3f, 0, H};\n' % (R, R))
    fid.write('Point(4) = {%.3f, %.3f, 0, H};\n\n' % (R, R))

    # 4.2) AIRFOILS' POINTS
    SUM = 0
    N[0] = copy.copy(
        SUM)  # N saves the number of points of the previous airfoils so that a continuous definition is possible.

    if nairfoils > 1:
        for p in range(1, nairfoils):
            SUM = SUM + len(x[p - 1])
            N[p] = copy.copy(SUM)

    for k in range(0, nairfoils):

        for k1 in range(0, len(x[k])):
            fid.write('Point(100 + %d) = {%.12f,%.12f,0,h}; \n' % (N[k] + k1, x[k][k1], y[k][k1]))

    # 4.3) ELLIPSE'S POINTS, LINES AND SURFACE
    for m in range(0, len(x_ellipse)):
        fid.write('Point(100000 + %d) = {%.12f,%.12f,0,h}; \n' % (m, x_ellipse[m], y_ellipse[m]))


    if len(gamma) == 2:
        fid.write(
            'Spline(5) = {(100000 + %d): (100000 + %d)};\n' % (gamma[0][0], gamma[1][0]))
        fid.write(
            'Spline(5 + 1) = {(100000 + %d):100359, 100000:(100000 + %d)};\n' % (gamma[1][0], gamma[0][0]))
        if ecc >= 0.01:
            fid.write('Transfinite Line{5} = %d Using Bump %.5f;\n' % (ellnodes[0], ecc / 4))
            fid.write('Transfinite Line{5 + 1} = %d Using Bump %.5f;\n' % (ellnodes[1], ecc / 4))
        else:
            fid.write('Transfinite Line{5} = %d Using Bump %.5f;\n' % (ellnodes[0], 0.0025))
            fid.write('Transfinite Line{5 + 1} = %d Using Bump %.5f;\n' % (ellnodes[1], 0.0025))
    else:
        for ell in range(0, len(gamma)):

            if ell < (len(gamma) - 1):

                if gamma[ell][1] == 'LE' or gamma[ell + 1][1] == 'LE':
                    fid.write(
                        'Spline(5 + %d) = {(100000 + %d): (100000 + %d)};\n' % (ell, gamma[ell][0], gamma[ell + 1][0]))

                    if ecc >= 0.01:
                        fid.write('Transfinite Line{5 + %d} = %d Using Bump %.5f;\n' % (ell, ellnodes[ell], ecc / 5))

                    else:
                        fid.write('Transfinite Line{5 + %d} = %d Using Bump %.5f;\n' % (ell, ellnodes[ell], 0.002))
                else:
                    fid.write(
                        'Spline(5 + %d) = {(100000 + %d): (100000 + %d)};\n' % (ell, gamma[ell][0], gamma[ell + 1][0]))
                    fid.write('Transfinite Line{5 + %d} = %d;\n' % (ell, ellnodes[ell]))

            else:

                if gamma[ell][1] == 'LE' or gamma[0][1] == 'LE':
                    fid.write('Spline(5 + %d) = {(100000 + %d):100359, 100000:(100000 + %d)};\n' % (
                        ell, gamma[ell][0], gamma[0][0]))

                    if ecc >= 0.01:
                        fid.write('Transfinite Line{5 + %d} = %d Using Bump %.5f;\n' % (ell, ellnodes[ell], ecc / 5))

                    else:
                        fid.write('Transfinite Line{5 + %d} = %d Using Bump %.5f;\n' % (ell, ellnodes[ell], 0.002))
                else:
                    fid.write('Spline(5 + %d) = {(100000 + %d):100359, 100000:(100000 + %d)};\n' % (
                        ell, gamma[ell][0], gamma[0][0]))
                    fid.write('Transfinite Line{5 + %d} = %d;\n' % (ell, ellnodes[ell]))

    fid.write('Curve Loop(2000) = {5, ')
    for j in range(1, len(gamma)):
        if j <= (len(gamma) - 2):
            fid.write('5 + %d, ' % j)
        else:
            fid.write('5 + %d};\n\n' % j)

    # 4.4) AIRFOILS' SPLINES, SURFACES AND STRUCTURED REGIONS
    for j in range(1, nairfoils + 1):
        i = j - 1

        fid.write('Spline(100*%d) = {100+%d:100+%d};\n' % (j, N[i], (N[i] + index_Below2[i])))
        fid.write(
            'Spline(100*%d + 1) = {100+%d:100+%d};\n' % (j, (N[i] + index_Below2[i]), (N[i] + index_Below1[i])))

        fid.write(
            'Spline(100*%d + 2) = {100+%d:100+%d};\n' % (j, (N[i] + index_Below1[i]), (N[i] + index_Above1[i])))

        fid.write(
            'Spline(100*%d + 3) = {100+%d:100+%d};\n' % (j, (N[i] + index_Above1[i]), (N[i] + index_Above2[i])))

        if TE[i] == 'Y':
            fid.write('Spline(100*%d + 4) = {100+%d:100+%d};\n' % (j, N[i] + index_Above2[i], N[i] + len(x[i]) - 1))
            fid.write('Line(100*%d + 5) = {100+%d, 100+%d};\n' % (j, N[i] + len(x[i]) - 1, N[i]))

        else:
            fid.write('Spline(100*%d + 4) = {100+%d:100+%d, 100+%d};\n' % (
                j, N[i] + index_Above2[i], N[i] + len(x[i]) - 1, N[i]))

        fid.write('Transfinite Line{100*%d} = %d Using Progression 1.015;\n' % (j, curveTE_L[i]))
        fid.write('Transfinite Line{100*%d + 1} = %d;\n' % (j, curveL[i]))

        fid.write('Transfinite Line{100*%d + 2} = %d Using Bump 5;\n' % (j, curveLE[i]))

        fid.write('Transfinite Line{100*%d + 3} = %d;\n' % (j, curveU[i]))
        fid.write('Transfinite Line{100*%d + 4} = %d Using Progression 0.985;\n' % (j, curveTE_U[i]))


        if TE[i] == 'Y':
            fid.write('Transfinite Line{100*%d + 5} = %d;\n' % (j, int(TE_len[i] * 20000 * nodes / crel[i])))
        else:
            pass


    # 4.5) FARFIELD
    fid.write('Circle(1) = {1,99,2};\n')
    fid.write('Line(2) = {2,3};\n')
    fid.write('Line(3) = {3,4};\n')
    fid.write('Line(4) = {4,1};\n')

    nodes_farfield = [int(30 / max(crel) * nodes), int(25 / max(crel) * nodes), int(50 / max(crel) * nodes)]
    fid.write('Transfinite Line{1} = %d Using Bump 2;\n' % nodes_farfield[0])
    fid.write('Transfinite Line{2,-4} = %d Using Progression 0.995;\n' % nodes_farfield[1])
    fid.write('Transfinite Line{3} = %d Using Bump 3.5;\n\n' % nodes_farfield[2])

    fid.write('Curve Loop(1000) = {1,2,3,4};\n')

    fid.write('Plane Surface(1000) = {1000, 2000}; \n\n')

    # 4.6) ELLIPSE SURFACE DEFINITION & CONCLUSION OF REFINING POINTS and AIRFOILS' SURFACES
    # Its definition was possible only after the definition of each airfoil curve loop.

    for j in range(1, nairfoils + 1):
        if TE[j - 1] == 'Y':
            if abs(x[j - 1][0] - x[j - 1][len(x[j - 1]) - 1]) <= 1e-4:
                m1[j - 1] = 0
            else:
                m1[j - 1] = - 1 / ((y[j - 1][0] - y[j - 1][len(y[j - 1]) - 1]) / (
                        x[j - 1][0] - x[j - 1][len(x[j - 1]) - 1]))
            x1[j - 1] = x[j - 1][0] + crel[j-1]/100 * mt.cos(mt.atan(m1[j - 1]))
            y1[j - 1] = y[j - 1][0] + crel[j-1]/100 * mt.sin(mt.atan(m1[j - 1]))
            x2[j - 1] = x[j - 1][len(x[j - 1]) - 1] + crel[j-1]/100 * mt.cos(mt.atan(m1[j - 1]))
            y2[j - 1] = y[j - 1][len(x[j - 1]) - 1] + crel[j-1]/100 * mt.sin(mt.atan(m1[j - 1]))
            # fid.write('Point(10*%d) = {%.12f,%.12f,0,h/4}; \n' % (j, x1[j - 1], y1[j - 1]))
            # fid.write('Point(10*%d + 1) = {%.12f,%.12f,0,h/4}; \n' % (j, x2[j - 1], y2[j - 1]))

            fid.write('Curve Loop(100*%d) = {100*%d, 100*%d + 1, 100*%d + 2, 100*%d + 3, 100*%d + 4, 100*%d + 5};\n' % (
                j, j, j, j, j, j, j))


        else:
            m1[j - 1] = [None] * 2

            m1[j - 1][0] = (y[j - 1][0] - y[j - 1][1]) / (x[j - 1][0] - x[j - 1][1])
            m1[j - 1][1] = (y[j - 1][len(x[j - 1]) - 1] - y[j - 1][len(x[j - 1]) - 2]) / (
                    x[j - 1][len(x[j - 1]) - 1] - x[j - 1][len(x[j - 1]) - 2])
            x1[j - 1] = x[j - 1][0] + crel[j-1]/100 * mt.cos(mt.atan(m1[j - 1][0]))
            y1[j - 1] = y[j - 1][0] + crel[j-1]/100 * mt.sin(mt.atan(m1[j - 1][0]))
            x2[j - 1] = x[j - 1][0] + crel[j-1]/100 * mt.cos(mt.atan(m1[j - 1][1]))
            y2[j - 1] = y[j - 1][0] + crel[j-1]/100 * mt.sin(mt.atan(m1[j - 1][1]))
            # fid.write('Point(10*%d) = {%.12f,%.12f,0,h/4}; \n' % (j, x1[j - 1], y1[j - 1]))
            # fid.write('Point(10*%d + 1) = {%.12f,%.12f,0,h/4}; \n' % (j, x2[j - 1], y2[j - 1]))
            fid.write(
                'Curve Loop(100*%d) = {100*%d, 100*%d + 1, 100*%d + 2, 100*%d + 3, 100*%d + 4};\n' % (j, j, j, j, j, j))

    fid.write('Plane Surface(2000) = {2000')
    for j in range(1, nairfoils + 1):
        if nairfoils == 1:
            fid.write(', 100}; \n\n')
        elif j != nairfoils:
            fid.write(', 100*%d' % j)
        else:
            fid.write(', 100*%d};\n\n' % j)

    # for i in range(0, nairfoils):
    #     if TE[i] == 'Y':
    #         pass
    #     else:
    #         fid.write('Point{10*%d} In Surface {2000}; \n' % (i + 1))
    #         fid.write('Point{10*%d + 1} In Surface {2000}; \n' % (i + 1))

    for j in range(0, len(gamma)):
        fid.write('Point{100000 + %d} In Surface {1000}; \n' % gamma[j][0])
        fid.write('Point{100000 + %d} In Surface {2000}; \n' % gamma[j][0])

    # 4.7) REFINING POINTS (WAKE)
    for k in range(0, nairfoils):

        if nairfoils == 1 or k == (nairfoils - 1):
            mindiff4 = 1000000

            if nairfoils == 1:

                if TE[k] == 'Y':
                    q2 = (y[k][0] + y[k][len(x) - 1]) / 2 - mt.tan(alpha[k] * mt.pi / 180) * (
                            x[k][0] + x[k][len(x) - 1]) / 2
                else:
                    q2 = y[k][0] - mt.tan(alpha[k] * mt.pi / 180) * x[k][0]

                for i in range(0, len(x_ellipse)):
                    if x_ellipse[i] >= x[k][0]:
                        diff4 = abs(y_ellipse[i] - (mt.tan(alpha[k] * mt.pi / 180) * x_ellipse[i] + q2))

                        if diff4 <= mindiff4:
                            mindiff4 = diff4
                            indicator[k] = i

            else:
                if TE[k] == 'Y':
                    q2 = (y[k][0] + y[k][len(x) - 1]) / 2 - mt.tan(
                        (- alpha[k] + alpha[ref_airfoil - 1]) * mt.pi / 180) * (x[k][0] + x[k][len(x) - 1]) / 2
                else:
                    q2 = y[k][0] - mt.tan((- alpha[k] + alpha[ref_airfoil - 1]) * mt.pi / 180) * x[k][0]

                for i in range(0, len(x_ellipse)):
                    if x_ellipse[i] >= x[k][0]:
                        diff4 = abs(y_ellipse[i] - (
                                    mt.tan((-alpha[k] + alpha[ref_airfoil - 1]) * mt.pi / 180) * x_ellipse[i] + q2))

                        if diff4 <= mindiff4:
                            mindiff4 = diff4
                            indicator[k] = i

            if TE[k] == 'Y':
                if alpha[k] >= 0:
                    x_TE1[k] = (x[k][0] + x[k][1] + x[k][len(x[k]) - 1] + x[k][len(x[k]) - 2] * 2 + x_ellipse[
                        indicator[k]]) / 6
                    y_TE1[k] = (y[k][0] + y[k][1] + y[k][len(x[k]) - 1] + y[k][len(x[k]) - 2] * 2 + y_ellipse[
                        indicator[k]]) / 6
                else:
                    x_TE1[k] = (x[k][0] + x[k][1] * 2 + x[k][len(x[k]) - 1] + x[k][len(x[k]) - 2] + x_ellipse[
                        indicator[k]]) / 6
                    y_TE1[k] = (y[k][0] + y[k][1] * 2 + y[k][len(x[k]) - 1] + y[k][len(x[k]) - 2] + y_ellipse[
                        indicator[k]]) / 6
                x_TE2[k] = (x[k][0] + x[k][1] + x[k][len(x[k]) - 1] + x[k][len(x[k]) - 2] + x_ellipse[indicator[k]]) / 5
                y_TE2[k] = (y[k][0] + y[k][1] + y[k][len(x[k]) - 1] + y[k][len(x[k]) - 2] + y_ellipse[indicator[k]]) / 5
                x_TE3[k] = (x[k][0] + x[k][len(x[k]) - 1] + x_ellipse[indicator[k]]) / 3
                y_TE3[k] = (y[k][0] + y[k][len(x[k]) - 1] + y_ellipse[indicator[k]]) / 3
                x_TE4[k] = ((x[k][0] + x[k][len(x[k]) - 1]) / 2 + x_ellipse[indicator[k]]) / 2
                y_TE4[k] = ((y[k][0] + y[k][len(x[k]) - 1]) / 2 + y_ellipse[indicator[k]]) / 2
                x_TE5[k] = ((x[k][0] + x[k][len(x[k]) - 1]) / 2 + x_ellipse[indicator[k]] * 2) / 3
                y_TE5[k] = ((y[k][0] + y[k][len(x[k]) - 1]) / 2 + y_ellipse[indicator[k]] * 2) / 3
                d[k] = mt.sqrt(
                    ((x[k][0] + x[k][len(x[k]) - 1]) / 2 - x_ellipse[indicator[k]]) ** 2 + (
                            (y[k][0] + y[k][len(x[k]) - 1]) / 2 - y_ellipse[indicator[k]]) ** 2)
                x_TE6[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] / 5
                y_TE6[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2
                x_TE7[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] / 3
                y_TE7[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2
                x_TE8[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] * 2 / 3
                y_TE8[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2

            else:
                x_TE1[k] = (x[k][0] * 5 + x_ellipse[indicator[k]]) / 6
                y_TE1[k] = (y[k][0] * 5 + y_ellipse[indicator[k]]) / 6
                x_TE2[k] = (x[k][0] * 4 + x_ellipse[indicator[k]]) / 5
                y_TE2[k] = (y[k][0] * 4 + y_ellipse[indicator[k]]) / 5
                x_TE3[k] = (x[k][0] * 2 + x_ellipse[indicator[k]]) / 3
                y_TE3[k] = (y[k][0] * 2 + y_ellipse[indicator[k]]) / 3
                x_TE4[k] = (x[k][0] + x_ellipse[indicator[k]]) / 2
                y_TE4[k] = (y[k][0] + y_ellipse[indicator[k]]) / 2
                x_TE5[k] = (x[k][0] + x_ellipse[indicator[k]] * 2) / 3
                y_TE5[k] = (y[k][0] + y_ellipse[indicator[k]] * 2) / 3
                d[k] = mt.sqrt(
                    (x[k][0] - x_ellipse[indw2]) ** 2 + (
                            y[k][0] - y_ellipse[indw2]) ** 2)
                x_TE6[k] = x[k][0] + d[k] / 5
                y_TE6[k] = y[k][0]
                x_TE7[k] = x[k][0] + d[k] / 3
                y_TE7[k] = y[k][0]
                x_TE8[k] = x[k][0] + d[k] * 2 / 3
                y_TE8[k] = y[k][0]

        else:
            MINDIFF = 10000
            if x[k + 1][index_sx[k + 1]] >= max(x[k][0], x[k][len(x[k]) - 1]):
                X1 = x[k][np.argmax(x[k])]
                Y1 = y[k][np.argmax(x[k])]
                # X2 = xc[k + 1][index_sx[k + 1]]
                # Y2 = yc[k + 1][index_sx[k + 1]]

            elif x[k + 1][index_sx[k + 1]] < max(x[k][0], x[k][len(x[k]) - 1]) and y[k + 1][index_sx[k + 1]] >= \
                    y[k][len(x[k]) - 1]:
                X1 = x[k][len(x[k]) - 1]
                Y1 = y[k][len(x[k]) - 1]
                # X2 = xc[k + 1][index_Below1[k + 1]]
                # Y2 = yc[k + 1][index_Below1[k + 1]]

            else:
                X1 = x[k][0]
                Y1 = y[k][0]
                # X2 = xc[k + 1][index_Above1[k + 1]]
                # Y2 = yc[k + 1][index_Above1[k + 1]]

            for i in range(0, len(x[k + 1])):
                DIFF = mt.sqrt((X1 - x[k + 1][i]) ** 2 + (Y1 - y[k + 1][i]) ** 2)
                if DIFF <= MINDIFF:
                    MINDIFF = DIFF
                    ind_gap = i

            x_TE1[k] = (X1 * 5 + x[k + 1][ind_gap]) / 6
            y_TE1[k] = (Y1 * 5 + y[k + 1][ind_gap]) / 6

            if ind_gap >= 2:
                x_TE2[k] = (X1 * 4 + x[k + 1][ind_gap - 1]) / 5
                y_TE2[k] = (Y1 * 4 + y[k + 1][ind_gap - 1]) / 5
                x_TE3[k] = (X1 * 2 + x[k + 1][ind_gap - 1]) / 3
                y_TE3[k] = (Y1 * 2 + y[k + 1][ind_gap - 1]) / 3
                x_TE4[k] = (X1 + x[k + 1][ind_gap - 2]) / 2
                y_TE4[k] = (Y1 + y[k + 1][ind_gap - 2]) / 2
                x_TE5[k] = (X1 + x[k + 1][ind_gap - 2] * 2) / 3
                y_TE5[k] = (Y1 + y[k + 1][ind_gap - 2] * 2) / 3
            elif ind_gap >= 1:
                x_TE2[k] = (X1 * 4 + x[k + 1][ind_gap - 1]) / 5
                y_TE2[k] = (Y1 * 4 + y[k + 1][ind_gap - 1]) / 5
                x_TE3[k] = (X1 * 2 + x[k + 1][ind_gap - 1]) / 3
                y_TE3[k] = (Y1 * 2 + y[k + 1][ind_gap - 1]) / 3
                x_TE4[k] = (X1 + x[k + 1][ind_gap - 1]) / 2
                y_TE4[k] = (Y1 + x[k + 1][ind_gap - 1]) / 2
                x_TE5[k] = (X1 + y[k + 1][ind_gap - 1] * 2) / 3
                y_TE5[k] = (Y1 + y[k + 1][ind_gap - 1] * 2) / 3
            else:
                x_TE2[k] = (X1 * 4 + x[k + 1][ind_gap]) / 5
                y_TE2[k] = (Y1 * 4 + y[k + 1][ind_gap]) / 5
                x_TE3[k] = (X1 * 2 + x[k + 1][ind_gap]) / 3
                y_TE3[k] = (Y1 * 2 + y[k + 1][ind_gap]) / 3
                x_TE4[k] = (X1 + x[k + 1][ind_gap]) / 2
                y_TE4[k] = (Y1 + y[k + 1][ind_gap]) / 2
                x_TE5[k] = (X1 + x[k + 1][ind_gap] * 2) / 3
                y_TE5[k] = (Y1 + y[k + 1][ind_gap] * 2) / 3

            if TE[k] == 'Y':
                d[k] = mt.sqrt(
                    ((x[k][0] + x[k][len(x[k]) - 1]) / 2 - x[k + 1][index_sx[k + 1]]) ** 2 + (
                            (y[k][0] + y[k][len(x[k]) - 1]) / 2 - y[k + 1][index_sx[k + 1]]) ** 2)
                x_TE6[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] / 5
                y_TE6[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2
                x_TE7[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] / 3
                y_TE7[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2
                x_TE8[k] = (x[k][0] + x[k][len(x[k]) - 1]) / 2 + d[k] * 2 / 3
                y_TE8[k] = (y[k][0] + y[k][len(x[k]) - 1]) / 2

            else:
                d[k] = mt.sqrt(
                    (x[k][0] - x[k + 1][index_sx[k + 1]]) ** 2 + (
                            y[k][0] - y[k + 1][index_sx[k + 1]]) ** 2)
                x_TE6[k] = x[k][0] + d[k] / 5
                y_TE6[k] = y[k][0]
                x_TE7[k] = x[k][0] + d[k] / 3
                y_TE7[k] = y[k][0]
                x_TE8[k] = x[k][0] + d[k] * 2 / 3
                y_TE8[k] = y[k][0]

        # fid.write('Point(10*%d + 2) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE1[k], y_TE1[k]))
        # fid.write('Point{10*%d + 2} In Surface {2000}; \n' % (k + 1))
        # fid.write('Point(10*%d + 3) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE2[k], y_TE2[k]))
        # fid.write('Point{10*%d + 3} In Surface {2000}; \n' % (k + 1))
        # fid.write('Point(10*%d + 4) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE3[k], y_TE3[k]))
        # fid.write('Point{10*%d + 4} In Surface {2000}; \n' % (k + 1))
        # fid.write('Point(10*%d + 5) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE4[k], y_TE4[k]))
        # fid.write('Point{10*%d + 5} In Surface {2000}; \n' % (k + 1))
        # fid.write('Point(10*%d + 6) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE5[k], y_TE5[k]))
        # fid.write('Point{10*%d + 6} In Surface {2000}; \n' % (k + 1))
        # fid.write('Point(10*%d + 7) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE6[k], y_TE6[k]))
        # fid.write('Point{10*%d + 7} In Surface {2000}; \n' % (k + 1))
        # fid.write('Point(10*%d + 8) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE7[k], y_TE7[k]))
        # fid.write('Point{10*%d + 8} In Surface {2000}; \n' % (k + 1))
        # fid.write('Point(10*%d + 9) = {%.12f,%.12f,0,h/2}; \n' % (k + 1, x_TE8[k], y_TE8[k]))
        # fid.write('Point{10*%d + 9} In Surface {2000}; \n' % (k + 1))


    # 4.8) PHYSICAL SURFACE AND LINES
    fid.write('Physical Surface(1) = {1000, 2000};\n\n')
    print(TE)
    for j in range(1, nairfoils + 1):
        if TE[j - 1] == 'N':
            fid.write('Physical Line("airfoil%d") = {100*%d, 100*%d + 1, 100*%d + 2, 100*%d + 3, 100*%d + 4};\n' % (
                j, j, j, j, j, j))
        else:
            fid.write(
                'Physical Line("airfoil%d") = {100*%d, 100*%d + 1, 100*%d + 2, 100*%d + 3, 100*%d + 4, 100*%d + 5};\n' % (
                    j, j, j, j, j, j, j))

    fid.write('Physical Line("farfield") = {1,2,3,4};\n')

    # File .su2 generation.
    fid.write('Mesh.Format = 42;\n')
    if Mesh_Algo == 6:
        fid.write('Mesh.Algorithm = 6;\n')
    elif Mesh_Algo == 'DEFAULT' or Mesh_Algo == 5:
        fid.write('Mesh.Algorithm = 5;\n')
    else:
        fid.write('Mesh.Algorithm = %d;\n' % Mesh_Algo)

    ll = 0
    if external_pts == 'YES':
        x_circle_ext, y_circle_ext = PointsInCircum(min(max(sel_data.b * semicircle_dimension, sel_data.a * semicircle_dimension), 200 * max(crel)), 360)
        i = 0
        while i < 360:
            fid.write('Point(100000 + %d) = {%.12f,%.12f,0,%f}; \n' % (len(x_ellipse) + ll, x_circle_ext[i], y_circle_ext[i], h * semicircle_elem_factor))
            fid.write('Point{100000 + %d} In Surface {1000}; \n' % (len(x_ellipse) + ll))
            ll = ll + 1
            i = i + 10

    else:
        pass

    fid.write('Mesh 2;\n')

    fid.close()

    return []
