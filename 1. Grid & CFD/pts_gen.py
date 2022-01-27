import numpy as np
import math as mt
import copy
from itertools import chain
import matplotlib.pyplot as plt


def airfoil_pts(x, y, nairfoils, index_Above1, index_Below1, index_Above2, index_Below2, crel, airfoil_nodes):

    index_ymin, index_ymax, index_sx, index_dx = [None] * nairfoils, [None] * nairfoils, [
        None] * nairfoils, [None] * nairfoils
    y_ymin, y_ymax, x_sx, x_dx = [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [None] * nairfoils
    curveU, curveL, curveLE, curveTE_L, curveTE_U, = [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [None] * nairfoils
    lengthLE, lengthUPPER, lengthLOWER, lengthTE_U, lengthTE_L = [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [None] * nairfoils

    for i in range(0, nairfoils):
        index_ymin[i] = np.argmin(y[i])
        index_ymax[i] = np.argmax(y[i])
        index_sx[i] = np.argmin(x[i])
        index_dx[i] = np.argmax(x[i])
        y_ymin[i] = np.min(y[i])
        y_ymax[i] = np.max(y[i])
        x_sx[i] = np.min(x[i])
        x_dx[i] = np.max(x[i])

        lengthTE_L[i] = sum([mt.sqrt((x[i][j1] - x[i][j1 + 1]) ** 2 + (y[i][j1] - y[i][j1 + 1]) ** 2) for j1 in
                          range(0, index_Below2[i])])
        lengthLOWER[i] = sum([mt.sqrt((x[i][j2] - x[i][j2 + 1]) ** 2 + (y[i][j2] - y[i][j2 + 1]) ** 2) for j2 in
                        range(index_Below2[i], index_Below1[i])])
        lengthLE[i] = sum([mt.sqrt((x[i][j3] - x[i][j3 + 1]) ** 2 + (y[i][j3] - y[i][j3 + 1]) ** 2) for j3 in
                        range(index_Below1[i], index_Above1[i])])
        lengthUPPER[i] = sum([mt.sqrt((x[i][j4] - x[i][j4 + 1]) ** 2 + (y[i][j4] - y[i][j4 + 1]) ** 2) for j4 in
                        range(index_Above1[i], index_Above2[i])])
        lengthTE_U[i] = sum([mt.sqrt((x[i][j5] - x[i][j5 + 1]) ** 2 + (y[i][j5] - y[i][j5 + 1]) ** 2) for j5 in
                          range(index_Above2[i], len(x[i]) - 1)])

        curveLE[i] = int(airfoil_nodes[i][0] * lengthLE[i] / crel[i])
        curveU[i] = int(airfoil_nodes[i][1] * lengthUPPER[i] / crel[i])
        curveL[i] = int(airfoil_nodes[i][2] * lengthLOWER[i] / crel[i])
        curveTE_U[i] = int(airfoil_nodes[i][3] * lengthTE_U[i] / crel[i])
        curveTE_L[i] = int(airfoil_nodes[i][4] * lengthTE_L[i] / crel[i])

    PDX = np.argmax(x_dx)
    PSX = np.argmin(x_sx)

    return y_ymin, y_ymax, x_sx, x_dx, index_ymin, index_ymax, index_sx, index_dx, curveLE, curveU, curveL, curveTE_U, curveTE_L, PSX, PDX, lengthLE, lengthUPPER, lengthLOWER, lengthTE_U, lengthTE_L


def grid_pts(x, y, nairfoils, BL_thickness, index_sx, index_Above1, index_Below1, index_Above2, index_Below2, TE):
    pointer = [0] * nairfoils
    xc, yc, theta = [None] * nairfoils, [None] * nairfoils, [None] * nairfoils
    lengthLE_C, lengthUPPER_C, lengthLOWER_C, lengthTE_U_C, lengthTE_L_C = [None] * nairfoils, [None] * nairfoils, [
        None] * nairfoils, [None] * nairfoils, [None] * nairfoils

    for t in range(0, nairfoils):

        xc[t], yc[t], theta[t] = [None] * len(x[t]), [None] * len(y[t]), [None] * len(x[t])

        for k in range(0, len(x[t])):

            if k < index_sx[t]:

                if x[t][k + 1] != x[t][k]:
                    X1 = copy.copy(x[t][k + 1])
                    Y1 = copy.copy(y[t][k + 1])
                    X2 = copy.copy(x[t][k])
                    Y2 = copy.copy(y[t][k])

                else:
                    j = k + 2
                    while x[t][j] == x[t][k]:
                        j = j + 1

                    X1 = copy.copy(x[t][j])
                    Y1 = copy.copy(y[t][j])
                    X2 = copy.copy(x[t][k])
                    Y2 = copy.copy(y[t][k])

                xc[t][k] = x[t][k] + BL_thickness[t] * mt.sin(mt.atan((Y2 - Y1) / (X2 - X1)))
                yc[t][k] = y[t][k] - BL_thickness[t] * mt.cos(mt.atan((Y2 - Y1) / (X2 - X1)))

            elif k == index_sx[t]:

                if x[t][k + 1] != x[t][k - 1]:
                    X1 = copy.copy(x[t][k + 1])
                    Y1 = copy.copy(y[t][k + 1])
                    X2 = copy.copy(x[t][k - 1])
                    Y2 = copy.copy(y[t][k - 1])
                else:
                    j = k + 2
                    while x[t][j] == x[t][k]:
                        j = j + 1

                    X1 = copy.copy(x[t][j])
                    Y1 = copy.copy(y[t][j])
                    X2 = copy.copy(x[t][k - 1])
                    Y2 = copy.copy(y[t][k - 1])

                xc[t][k] = x[t][k] - abs(BL_thickness[t] * mt.sin(mt.atan((Y2 - Y1) / (X2 - X1))))
                yc[t][k] = y[t][k] - abs(BL_thickness[t] * mt.cos(mt.atan((Y2 - Y1) / (X2 - X1))))

            else:

                if x[t][k - 1] != x[t][k]:
                    X1 = copy.copy(x[t][k - 1])
                    Y1 = copy.copy(y[t][k - 1])
                    X2 = copy.copy(x[t][k])
                    Y2 = copy.copy(y[t][k])

                else:
                    j = k - 2
                    while x[t][j] == x[t][k]:
                        j = j - 1

                    X1 = copy.copy(x[t][j])
                    Y1 = copy.copy(y[t][j])
                    X2 = copy.copy(x[t][k])
                    Y2 = copy.copy(y[t][k])

                xc[t][k] = x[t][k] - BL_thickness[t] * mt.sin(mt.atan((Y2 - Y1) / (X2 - X1)))
                yc[t][k] = y[t][k] + BL_thickness[t] * mt.cos(mt.atan((Y2 - Y1) / (X2 - X1)))

        if TE[t] == 'N':
            if x[t][len(x[t]) - 1] != x[t][0]:
                X1 = copy.copy(x[t][len(x[t]) - 1])
                Y1 = copy.copy(y[t][len(x[t]) - 1])
                X2 = copy.copy(x[t][0])
                Y2 = copy.copy(y[t][0])

            else:
                j = len(x[t]) - 2
                while x[t][j] == x[t][0]:
                    j = j - 1

                X1 = copy.copy(x[t][j])
                Y1 = copy.copy(y[t][j])
                X2 = copy.copy(x[t][0])
                Y2 = copy.copy(y[t][0])

            xc[t].append(x[t][0] - BL_thickness[t] * mt.sin(mt.atan((Y2 - Y1) / (X2 - X1))))
            yc[t].append(y[t][0] + BL_thickness[t] * mt.cos(mt.atan((Y2 - Y1) / (X2 - X1))))
            pointer[t] = 1
        else:
            pass

        lengthTE_L_C[t] = sum([mt.sqrt((xc[t][j1] - xc[t][j1 + 1]) ** 2 + (yc[t][j1] - yc[t][j1 + 1]) ** 2) for j1 in
                               range(0, index_Below2[t])])
        lengthLOWER_C[t] = sum([mt.sqrt((xc[t][j2] - xc[t][j2 + 1]) ** 2 + (yc[t][j2] - yc[t][j2 + 1]) ** 2) for j2 in
                                range(index_Below2[t], index_Below1[t])])
        lengthLE_C[t] = sum([mt.sqrt((xc[t][j3] - xc[t][j3 + 1]) ** 2 + (yc[t][j3] - yc[t][j3 + 1]) ** 2) for j3 in
                             range(index_Below1[t], index_Above1[t])])
        lengthUPPER_C[t] = sum([mt.sqrt((xc[t][j4] - xc[t][j4 + 1]) ** 2 + (yc[t][j4] - yc[t][j4 + 1]) ** 2) for j4 in
                                range(index_Above1[t], index_Above2[t])])
        lengthTE_U_C[t] = sum([mt.sqrt((xc[t][j5] - xc[t][j5 + 1]) ** 2 + (yc[t][j5] - yc[t][j5 + 1]) ** 2) for j5 in
                               range(index_Above2[t], len(x[t]) - 1)])

    return xc, yc, theta, lengthLE_C, lengthUPPER_C, lengthLOWER_C, lengthTE_U_C, lengthTE_L_C, pointer


class ell:
    def __init__(self, a, b, area, enclose):
        self.a = a
        self.b = b
        self.area = area
        self.enclose = enclose


def ext_pts_VISCOUS(x, y, xc, yc, nairfoils, index_sx, PSX, PDX, crel, ellipse_dimension):

    def plot_ellipse_axis(x, y, b):
        # plotting the actual points as scatter plot
        plt.scatter(x, y, color="m",
                    marker="o", s=30)

        # predicted response vector
        y_pred = b[0] + b[1] * x

        # plotting the regression line
        plt.plot(x, y_pred, color="g")

        # putting labels
        plt.xlabel('x')
        plt.ylabel('y')

        # function to show plot
        plt.savefig('airfoil_image2.png')
        plt.close()

    x_1D = list(chain.from_iterable([list(chain.from_iterable(x)), list(chain.from_iterable(xc))]))
    y_1D = list(chain.from_iterable([list(chain.from_iterable(y)), list(chain.from_iterable(yc))]))
    mean_x = np.mean(np.asarray(x_1D))
    mean_y = np.mean(np.asarray(y_1D))
    sum_x_centre = 0
    sum_y_centre = 0

    for j in range(0, nairfoils):
        sum_x_centre = sum_x_centre + (x[j][0] + x[j][index_sx[j]]) / 2 * crel[j]
        sum_y_centre = sum_y_centre + (y[j][0] + y[j][index_sx[j]]) / 2 * crel[j]

    centre_x = sum_x_centre / sum(crel)
    centre_y = sum_y_centre / sum(crel)

    # estimating coefficients
    b = [y[PDX][0] - (y[PDX][0] - y[PSX][index_sx[PSX]])/(x[PDX][0] - x[PSX][index_sx[PSX]]) * x[PDX][0], (y[PDX][0] - y[PSX][index_sx[PSX]])/(x[PDX][0] - x[PSX][index_sx[PSX]])]
    # plotting regression line
    #plot_ellipse_axis(np.array(x_1D), np.array(y_1D), b)
    tau = mt.atan(b[1])

    dist1 = np.max(abs(centre_x - x_1D))
    dist2 = np.max(abs(centre_y - y_1D))

    a = np.linspace(dist1, 5*dist1, num=50)
    b = np.linspace(dist2, 5*dist2, num=50)

    x_ellipse, y_ellipse, dat = [None] * len(a) * len(b), [None] * len(a) * len(b), [None] * len(a) * len(b)

    ii = 0
    N = 360

    for i in range(0, len(a)):
        for k in range(0, len(b)):

            quest = 1

            for j in range(len(x_1D)):

                xc = x_1D[j] - centre_x
                yc = y_1D[j] - centre_y

                xct = xc * np.cos(tau) + yc * np.sin(tau)
                yct = xc * np.sin(tau) - yc * np.cos(tau)

                rad_cc = (xct ** 2 / (a[i]) ** 2) + (yct ** 2 / b[k] ** 2)

                if rad_cc < ellipse_dimension:
                    pass
                else:
                    quest = 0

            dat[ii] = ell(a[i], b[k], mt.pi * a[i] * b[i], quest)

            ii = ii + 1

    minA = 100000000000
    for ii2 in range(0, len(dat)):
        if mean_x >= mean_y:
            if dat[ii2].area < minA and dat[ii2].enclose == 1 and dat[ii2].a >= dat[ii2].b:
                minA = dat[ii2].area
                index = ii2
        else:
            if dat[ii2].area < minA and dat[ii2].enclose == 1 and dat[ii2].a < dat[ii2].b:
                minA = dat[ii2].area
                index = ii2

    sel_data = copy.copy(dat[index])

    t = np.linspace(0, N - 1, num=N)
    x_ellipse = [None] * N
    y_ellipse = [None] * N

    for jj in range(0, 360):
        x_ellipse[jj] = centre_x + mt.cos(tau) * sel_data.a * mt.cos(t[jj] * mt.pi / 180) - mt.sin(
            tau) * sel_data.b * mt.sin(t[jj] * mt.pi / 180)
        y_ellipse[jj] = centre_y + mt.sin(tau) * sel_data.a * mt.cos(t[jj] * mt.pi / 180) + mt.cos(
            tau) * sel_data.b * mt.sin(t[jj] * mt.pi / 180)

    if sel_data.a >= sel_data.b:
        ecc = mt.sqrt(sel_data.a ** 2 - sel_data.b ** 2) / sel_data.a

    else:
        ecc = mt.sqrt(sel_data.b ** 2 - sel_data.a ** 2) / sel_data.b

    return centre_x, centre_y, x_ellipse, y_ellipse, sel_data, tau, ecc


def ext_pts_EULER(x, y, nairfoils, index_sx, PSX, PDX, crel, ellipse_dimension):

    def plot_ellipse_axis(x, y, b):
        # plotting the actual points as scatter plot
        plt.scatter(x, y, color="m",
                    marker="o", s=30)

        # predicted response vector
        y_pred = b[0] + b[1] * x

        # plotting the regression line
        plt.plot(x, y_pred, color="g")

        # putting labels
        plt.xlabel('x')
        plt.ylabel('y')

        # function to show plot
        plt.savefig('airfoil_image2.png')
        plt.close()

    x_1D = list(chain.from_iterable(x))
    y_1D = list(chain.from_iterable(y))
    mean_x = np.mean(np.asarray(x_1D))
    mean_y = np.mean(np.asarray(y_1D))

    sum_x_centre = 0
    sum_y_centre = 0

    for j in range(0, nairfoils):
        sum_x_centre = sum_x_centre + (x[j][0] + x[j][index_sx[j]]) / 2 * crel[j]
        sum_y_centre = sum_y_centre + (y[j][0] + y[j][index_sx[j]]) / 2 * crel[j]

    centre_x = sum_x_centre / sum(crel)
    centre_y = sum_y_centre / sum(crel)

    # estimating coefficients
    b = [y[PDX][0] - (y[PDX][0] - y[PSX][index_sx[PSX]])/(x[PDX][0] - x[PSX][index_sx[PSX]]) * x[PDX][0], (y[PDX][0] - y[PSX][index_sx[PSX]])/(x[PDX][0] - x[PSX][index_sx[PSX]])]
    # plotting regression line
    #plot_ellipse_axis(np.array(x_1D), np.array(y_1D), b)
    tau = mt.atan(b[1])

    dist1 = np.max(abs(centre_x - x_1D))
    dist2 = np.max(abs(centre_y - y_1D))

    a = np.linspace(dist1, 5*dist1, num=50)
    b = np.linspace(dist2, 5*dist2, num=50)

    x_ellipse, y_ellipse, dat = [None] * len(a) * len(b), [None] * len(a) * len(b), [None] * len(a) * len(b)

    ii = 0
    N = 360

    for i in range(0, len(a)):
        for k in range(0, len(b)):

            quest = 1

            for j in range(len(x_1D)):

                xc = x_1D[j] - centre_x
                yc = y_1D[j] - centre_y

                xct = xc * np.cos(tau) + yc * np.sin(tau)
                yct = xc * np.sin(tau) - yc * np.cos(tau)

                rad_cc = (xct ** 2 / (a[i]) ** 2) + (yct ** 2 / b[k] ** 2)

                if rad_cc < ellipse_dimension:
                    pass
                else:
                    quest = 0

            dat[ii] = ell(a[i], b[k], mt.pi * a[i] * b[i], quest)

            ii = ii + 1

    minA = 100000000000
    for ii2 in range(0, len(dat)):
        if mean_x >= mean_y:
            if dat[ii2].area < minA and dat[ii2].enclose == 1 and dat[ii2].a >= dat[ii2].b:
                minA = dat[ii2].area
                index = ii2
        else:
            if dat[ii2].area < minA and dat[ii2].enclose == 1 and dat[ii2].a < dat[ii2].b:
                minA = dat[ii2].area
                index = ii2

    sel_data = copy.copy(dat[index])

    t = np.linspace(0, N - 1, num=N)
    x_ellipse = [None] * N
    y_ellipse = [None] * N

    for jj in range(0, 360):
        x_ellipse[jj] = centre_x + mt.cos(tau) * sel_data.a * mt.cos(t[jj] * mt.pi / 180) - mt.sin(
            tau) * sel_data.b * mt.sin(t[jj] * mt.pi / 180)
        y_ellipse[jj] = centre_y + mt.sin(tau) * sel_data.a * mt.cos(t[jj] * mt.pi / 180) + mt.cos(
            tau) * sel_data.b * mt.sin(t[jj] * mt.pi / 180)

    if sel_data.a >= sel_data.b:
        ecc = mt.sqrt(sel_data.a ** 2 - sel_data.b ** 2) / sel_data.a

    else:
        ecc = mt.sqrt(sel_data.b ** 2 - sel_data.a ** 2) / sel_data.b

    return centre_x, centre_y, x_ellipse, y_ellipse, sel_data, tau, ecc


