# IGP Parameterization Profiles
from getThickParamIGP import getThickParamIGP
import matplotlib.pyplot as plt
from camberline import camberlineIGP
from camberline import camberlineNACA4
from camberline import camberlineNACA5_REFL
from camberline import camberlineNACA5_NONREFL
from camberline import camberlineNACA4_MOD
import numpy as np
import math as mt
import sys


def getNpt(n_points, relChord):
    min_npt = 0.2 * n_points
    fctr = 4 / 3
    
    npt = fctr * relChord * n_points
    npt = round(npt)
    npt = min(npt, n_points)
    npt = max(npt, min_npt)

    return npt


def Finite_TE_correction(x, y):

    m_u = (y[len(x) - 1] - y[len(x) - 2]) / (x[len(x) - 1] - x[len(x) - 2])
    q_u = y[len(x) - 1] - m_u * x[len(x) - 1]
    m_l = (y[0] - y[1]) / (x[0] - x[1])
    q_l = y[0] - m_l * x[0]
    x_new = (q_l - q_u) / (m_u - m_l)
    y_new = m_u * x_new + q_u
    xl = x.tolist()
    xl.insert(0, x_new)
    yl = y.tolist()
    yl.insert(0, y_new)

    x = np.asarray(xl)
    y = np.asarray(yl)

    c = np.max(x) - np.min(x)
    x = x / c
    y = y / c

    x = x.tolist()
    y = y.tolist()
    for i in range(0, 5):
        x.pop(2)
        y.pop(2)

        x.pop(len(x) - 3)
        y.pop(len(y) - 3)


    x = np.asarray(x)
    y = np.asarray(y)

    return x, y


def airfoilIGP(param, n, TE):
    c1 = param[0]
    c2 = param[1]
    c3 = param[2]
    c4 = param[3]
    xt = param[4]
    T = param[5]
    rho0_bar = param[6]
    beta_te_bar = param[7]

    t1, t2, t3, t4, t5 = getThickParamIGP(xt, T, rho0_bar, beta_te_bar)

    xu, yu, xl, yl = np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n])

    for i in range(1, n + 1):
        k = 1 - 0.5 * (1 + np.cos(((i - 1) * np.pi) / (n - 1)))

        xc, yc = camberlineIGP(c1, c2, c3, c4, k)

        thickfun = lambda x: t1 * np.sqrt(x) + t2 * x + t3 * x ** 2 + t4 * x ** 3 + t5 * x ** 4

        t = thickfun(k)

        xu[i - 1] = xc
        yu[i - 1] = yc + 0.5 * t

        xl[i - 1] = xc
        yl[i - 1] = yc - 0.5 * t
    xu = xu.tolist()
    xu.pop(len(xu) - 1)
    yu = yu.tolist()
    yu.pop(len(yu) - 1)
    xu = np.asarray(xu)
    yu = np.asarray(yu)
    x = np.concatenate([np.flip(xl, 0), xu[1:]])
    y = np.concatenate([np.flip(yl, 0), yu[1:]])
    if TE == 'Y':
        print('IGP airfoils have been implemented with finite angle trailing edge only. As consequence, TE = "N" is imposed for each IGP airfoil.\n')
    TE = 'N'

    return x, y, TE


def airfoilNACA4(param, n, TE):
    m = param[0] / 100
    p = param[1] / 10
    t = param[2] / 100
    c = 1

    xu, yu, xl, yl = np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n])

    xc = [0] * n

    # if TE == 'N':
    #     xc2 = np.linspace(0, mt.pi / 2, n)
    #     for j in range(1, n):
    #         xc[j] = abs(mt.sin(xc2[j]))
    #     xc.reverse()
    #else:
    xc2 = np.linspace(mt.pi / 2, mt.pi, n)
    for j in range(1, n):
        xc[j] = pow(mt.cos(xc2[j]), 2)

    xc = np.asarray(xc)

    for i in range(0, n):
        yc, theta, yt = camberlineNACA4(xc[i], m, p, t, c)

        xu[i] = xc[i] - yt * mt.sin(theta)
        yu[i] = yc + yt * mt.cos(theta)
        xl[i] = xc[i] + yt * mt.sin(theta)
        yl[i] = yc - yt * mt.cos(theta)

    x = np.concatenate([np.flip(xl, 0), xu[1:]])
    y = np.concatenate([np.flip(yl, 0), yu[1:]])

    if TE == 'N':
        x, y = Finite_TE_correction(x, y)
    else:
        pass

    return x, y


def airfoilNACA5(param, n, TE):
    c = 1

    xu, yu, xl, yl = np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n])

    xc2 = np.linspace(mt.pi / 2, mt.pi, n)
    xc = [0] * n

    for j in range(1, n):
        xc[j] = pow(mt.cos(xc2[j]), 2)

    xc = np.asarray(xc)

    for i in range(0, n):
        if param[2] == 0:
            yc, theta, yt = camberlineNACA5_NONREFL(xc[i], c, param)

        elif param[2] == 1:
            yc, theta, yt = camberlineNACA5_REFL(xc[i], c, param)

        else:
            print("Error: accepted input about third parameter are only 0, 1 values. Please Check input parameters")
            sys.exit()

        xu[i] = xc[i] - yt * mt.sin(theta)
        yu[i] = yc + yt * mt.cos(theta)
        xl[i] = xc[i] + yt * mt.sin(theta)
        yl[i] = yc - yt * mt.cos(theta)

        x = np.concatenate([np.flip(xl, 0), xu[1:]])
        y = np.concatenate([np.flip(yl, 0), yu[1:]])

    if TE == 'N':
        x, y = Finite_TE_correction(x, y)
    else:
        pass

    return x, y


def airfoilNACA4_MOD(param, n, TE):
    c = 1
    m = param[0] / 100
    p = param[1] / 10
    t = param[2] / 100
    I = param[3]
    T = param[4] / 10

    xu, yu, xl, yl = np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n])

    xc2 = np.linspace(mt.pi / 2, mt.pi, n)
    xc = [0] * n

    for j in range(1, n):
        xc[j] = pow(mt.cos(xc2[j]), 2)

    xc = np.asarray(xc)

    for i in range(0, n):
        yc, theta, yt = camberlineNACA4_MOD(xc[i], c, m, p, t, I, T)
        xu[i] = xc[i] - yt * mt.sin(theta)
        yu[i] = yc + yt * mt.cos(theta)
        xl[i] = xc[i] + yt * mt.sin(theta)
        yl[i] = yc - yt * mt.cos(theta)

        x = np.concatenate([np.flip(xl, 0), xu[1:]])
        y = np.concatenate([np.flip(yl, 0), yu[1:]])

    if TE == 'N':
        x, y = Finite_TE_correction(x, y)
    else:
        pass

    return x, y


def airfoilNACA16(param, n, TE):
    c = 1
    m = param[0] / 100
    p = param[1] / 10
    t = param[2] / 100
    I = 4
    T = 5 / 10

    xu, yu, xl, yl = np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n])

    xc2 = np.linspace(mt.pi / 2, mt.pi, n)
    xc = [0] * n

    for j in range(1, n):
        xc[j] = pow(mt.cos(xc2[j]), 2)

    xc = np.asarray(xc)

    for i in range(0, n):
        yc, theta, yt = camberlineNACA4_MOD(xc[i], c, m, p, t, I, T)
        xu[i] = xc[i] - yt * mt.sin(theta)
        yu[i] = yc + yt * mt.cos(theta)
        xl[i] = xc[i] + yt * mt.sin(theta)
        yl[i] = yc - yt * mt.cos(theta)

        x = np.concatenate([np.flip(xl, 0), xu[1:]])
        y = np.concatenate([np.flip(yl, 0), yu[1:]])

    if TE == 'N':
        x, y = Finite_TE_correction(x, y)
    else:
        pass

    return x, y


def airfoilBICONVEX(param, n, TE):
    c = 1

    xu, yu, xl, yl = np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n])

    xc2 = np.linspace(mt.pi / 2, mt.pi, n)
    xc = [0] * n

    for j in range(1, n):
        xc[j] = pow(mt.cos(xc2[j]), 2)

    xc = np.asarray(xc)

    for i in range(0, n):
        xu[i] = xc[i]
        yu[i] = param[0] * (xu[i] - (xu[i] ** param[1]))
        xl[i] = xc[i]
        yl[i] = - yu[i]

        x = np.concatenate([np.flip(xl, 0), xu[1:]])
        y = np.concatenate([np.flip(yl, 0), yu[1:]])

    if TE == 'N':
        x, y = Finite_TE_correction(x, y)
    else:
        pass

    return x, y


def airfoilMY_FILE(i, TE):

    with open("MY_FILE_%d.dat" % i) as f:
        data = [[float(row.split()[0]), float(row.split()[1])] for row in f]
        data1 = [data[j] for j in range(0, int(len(data) / 2))]
        data2 = [data[j] for j in range(len(data)-1, int(len(data)/2), -1)]

        if (sum([data1[i][1] for i in range(0, len(data1))]) / len(data1)) >= (sum([data2[i][1] for i in range(0, len(data2))]) / len(data2)):
            data = data[::-1]
        else:
            pass

    x = np.asarray([data[i][0] for i in range(0, len(data))])
    y = np.asarray([data[i][1] for i in range(0, len(data))])

    if TE == 'N' and x[len(x) - 1] == x[0] and y[len(y) - 1] == y[0]:
        x = np.delete(x, len(x) - 1)
        y = np.delete(y, len(y) - 1)

    x = np.asarray(x)
    y = np.asarray(y)

    return x, y, TE


def airfoilBENZING(param, n, TE, i):
    print('For Benzing airfoils there is no available option for trailing edge. These airfoils type are printed with finite angle on trailing edge as default.\n')

    xu, yu, xl, yl, yc, yt, theta = np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n]), np.zeros([n])

    xc2 = np.linspace(mt.pi / 2, mt.pi, n)
    xc = [0] * n

    for j in range(1, n):
        xc[j] = (pow(mt.cos(xc2[j]), 2)) * 100

    xc = np.asarray(xc)

    if param[1] == 46:
        if param[0] == 92:
            S = 80
        elif param[0] == 122:
            S = 105
        elif param[0] == 152:
            S = 131

        F = 0.25
        B0 = 0.1975
        B1 = 0.559
        B2 = -2.876/1000
        B3 = -2.71/100000
        B4 = 0.00000001/10000000

        for i in range(0, n):
            yt[i] = S * (0.26 * ((xc[i] / 100) ** 0.5) - 0.22 * (xc[i] / 100) - 0.3516 * ((xc[i] / 100) ** 2) + 0.5 * (
                        (xc[i] / 100) ** 3) - 0.1885 * ((xc[i] / 100) ** 4))
            yc[i] = -F * (B0 + B1 * xc[i] + B2 * (xc[i] ** 2) + B3 * (xc[i] ** 3) + B4 * (xc[i] ** 4))
            theta[i] = mt.atan(B1 + 2 * B2 * xc[i] + 3 * B3 * (xc[i] ** 2) + 4 * B4 * (xc[i] ** 3))

            xl[i] = (xc[i] - yt[i] * mt.sin(theta[i])) / 100
            yl[i] = (yc[i] - yt[i] * mt.cos(theta[i])) / 100
            xu[i] = (xc[i] + yt[i] * mt.sin(theta[i])) / 100
            yu[i] = (yc[i] + yt[i] * mt.cos(theta[i])) / 100


    elif param[1] == 75:
        if param[0] == 92:
            S = 75
        elif param[0] == 122:
            S = 92
        elif param[0] == 152:
            S = 122

        F = 0.72
        B0 = 0.2101
        B1 = 0.33
        B2 = - 1.945 / 1000
        B3 = - 2.849 / 100000
        B4 = 1.47401 / 10000000

        for i in range(0, n):
            yt[i] = S * (0.26 * ((xc[i] / 100) ** 0.5) - 0.22 * (xc[i] / 100) - 0.3516 * ((xc[i] / 100) ** 2) + 0.5 * ((xc[i] / 100) ** 3) - 0.1885 * ((xc[i] / 100) ** 4))
            yc[i] = -F * (B0 + B1 * xc[i] + B2 * (xc[i] ** 2) + B3 * (xc[i] ** 3) + B4 * (xc[i] ** 4))
            theta[i] = mt.atan(B1 + 2*B2*xc[i] + 3*B3*(xc[i]**2) + 4*B4*(xc[i]**3))

            xu[i] = (xc[i] + yt[i] * mt.sin(theta[i])) / 100

            xl[i] = (xc[i] - yt[i] * mt.sin(theta[i])) / 100
            yl[i] = (yc[i] - yt[i] * mt.cos(theta[i])) / 100
            xu[i] = (xc[i] + yt[i] * mt.sin(theta[i])) / 100
            yu[i] = (yc[i] + yt[i] * mt.cos(theta[i])) / 100


    elif param[1] == 56:
        if param[0] == 94:
            s = 1
            X0 = 0.7632/100
            Y0 = 0.0808/100
            R0 = 0.8040/100
        elif param[0] == 124:
            s = 1.255
            X0 = 1.0833/100
            Y0 = 0.0971/100
            R0 = 1.1044/100
        elif param[0] == 154:
            s = 1.56
            X0 = 1.4295/100
            Y0 = -0.1508/100
            R0 = 1.4548/100

        F = 0.27
        A0 = 0.6952458
        A1 = 0.2640887
        A2 = -4.749312
        A3 = 2.277717
        A4 = -2.547127

        for i in range(0, n):
            yt[i] = s * (A0 + A1*xc[i] + A2/1000*(xc[i]**2) + A3/100000*(xc[i]**3) + A4/100000000 * (xc[i]**4))
            yc[i] = (0.2508 - 0.3300893 * xc[i] - 0.001945241 * (xc[i]**2) + 0.02849299/1000 * (xc[i]**3) + 0.00237/10000 * (xc[i]**4)) * F
            theta[i] = mt.atan(-0.3300893 - 0.001945241 * 2 * xc[i] + 3 * 0.02849299 / 1000 * (xc[i] ** 2) + 4 * (0.00237 / 10000) * (xc[i]**3))

            xl[i] = (xc[i] - yt[i] * mt.sin(theta[i])) / 100
            yl[i] = (yc[i] - yt[i] * mt.cos(theta[i])) / 100
            xu[i] = (xc[i] + yt[i] * mt.sin(theta[i])) / 100
            yu[i] = (yc[i] + yt[i] * mt.cos(theta[i])) / 100

        xu = xu.tolist()
        yu = yu.tolist()
        xl = xl.tolist()
        yl = yl.tolist()
        i = 0
        while i < len(xl):
            if xl[i] <= X0 + R0:
                xl.pop(i)
                yl.pop(i)
            else:
                i = i + 1
        j = 0
        while i < len(xu):
            if xu[j] <= X0 + R0:
                xu.pop(j)
                yu.pop(j)
            else:
                i = i + 1

        angle = np.linspace(mt.pi /2, mt.pi * 3/2, 30)
        sel_angle1, sel_angle2 = [], []

        for j in range(1, len(angle)):

            if (X0 + R0 * mt.cos(angle[j])) < xu[0] and R0 * mt.sin(angle[j]) >= 0:
                sel_angle1.append(angle[j])
            if (X0 + R0 * mt.cos(angle[j])) < xl[0] and R0 * mt.sin(angle[j]) <= 0:
                sel_angle2.append(angle[j])
        print(sel_angle1[1] * 180 / mt.pi, sel_angle2[len(sel_angle1) - 2] * 180 / mt.pi)
        sel_angle1.pop(1)
        sel_angle1.pop(2)
        sel_angle2.pop(len(sel_angle1) - 2)
        sel_angle2.pop(len(sel_angle1) - 3)

        sel_angle1.reverse()
        ref_u = [[xu[0], xu[1]], [yu[0], yu[1]]]
        ref_l = [[xl[0], xl[1]], [yl[0], yl[1]]]

        for k1 in range(0, len(sel_angle1)):
            if (ref_u[1][1] - ref_u[1][0]) / (ref_u[0][1] - ref_u[0][0]) <= ((Y0 + R0 * mt.sin(sel_angle1[k1])) - (Y0 + R0 * mt.sin(sel_angle1[k1 - 1]))) / ((X0 + R0 * mt.cos(sel_angle1[k1])) - (X0 + R0 * mt.cos(sel_angle1[k1 - 1]))):
                xu.insert(k1, (X0 + R0 * mt.cos(sel_angle1[k1])))
                yu.insert(k1, (Y0 + R0 * mt.sin(sel_angle1[k1])))

        for k2 in range(0, len(sel_angle2)):
            if (ref_l[1][1] - ref_l[1][0]) / (ref_l[0][1] - ref_l[0][0]) >= ((Y0 + R0 * mt.sin(sel_angle2[k2])) - (Y0 + R0 * mt.sin(sel_angle2[k2 - 1]))) / ((X0 + R0 * mt.cos(sel_angle2[k2])) - (X0 + R0 * mt.cos(sel_angle2[k2 - 1]))):
                xl.insert(k2, (X0 + R0 * mt.cos(sel_angle2[k2])))
                yl.insert(k2, (Y0 + R0 * mt.sin(sel_angle2[k2])))

        xl = np.asarray(xl)
        xu = np.asarray(xu)
        yl = np.asarray(yl)
        yu = np.asarray(yu)

    elif param[1] == 77:
        if len(str(param[0])) == 3:
            num1 = str(param[0] - int(str(param[0])[2]))[0]
            num2 = str(param[0] - int(str(param[0])[2]))[1]
            Sr = int(num1 + num2)
        elif len(str(param[0])) == 2:
            Sr = int(str(param[0] - int(str(param[0])[1]))[0])
        else:
            print('Error: the selected Benzing airfoil does not exist into the database. Please check airfoil number ('+ i +') parameters. \n')
            sys.exit()

        for i in range(0, n):
            yc[i] = -0.09717 - 0.33406 * xc[i] + 0.00955 * (xc[i]**2) - 0.00016 * (xc[i]**3) + 0.000001 * (xc[i]**4)
            yt[i] = 100 * ((Sr/20) * (0.2969 * ((xc[i]/100) ** 0.5) - 0.126 * (xc[i]/100) - 0.3516 * ((xc[i]/100) ** 2) + 0.2843 * ((xc[i]/100)**3) - 0.1015 * ((xc[i]/100) **4)))
            theta[i] = mt.atan(0.33406 - 2*0.00955*xc[i] + 3*0.00016*(xc[i]**2) - 4 * 0.000001*(xc[i]**3))

            xl[i] = (xc[i] - yt[i] * mt.sin(theta[i])) / 100
            yl[i] = (yc[i] - yt[i] * mt.cos(theta[i])) / 100
            xu[i] = (xc[i] + yt[i] * mt.sin(theta[i])) / 100
            yu[i] = (yc[i] + yt[i] * mt.cos(theta[i])) / 100

    elif param[1] == 115:  # TO BE FIXED! PROVIDED RELATIONS BY BENZING DO NOT WORK!
        if param[0] == 103:
            K = 1.55
        elif param[0] == 123:
            K = 1.75
        elif param[0] == 153:
            K = 2.25

        for i in range(0, n):
            yt[i] = K * (0.398 + 0.2299 * xc[i] - (5.242/1000) * (xc[i]**2) + (3.321/100000) * (xc[i]**3) + (4.09/100000000) * (xc[i]**4))
            yc[i] = 0.1 - 0.5406 * xc[i] + 0.009755 * (xc[i]**2) - 0.0000916 * (xc[i]**3) + 0.00000048011 * (xc[i]**4)
            theta[i] = mt.atan(0.5406 - 2*0.009755*xc[i] + 3*0.0000916*(xc[i]**2) - 4*0.0000004811*(xc[i]**3))

            xl[i] = (xc[i] - yt[i] * mt.sin(theta[i])) / 100
            yl[i] = (yc[i] - yt[i] * mt.cos(theta[i])) / 100
            xu[i] = (xc[i] + yt[i] * mt.sin(theta[i])) / 100
            yu[i] = (yc[i] + yt[i] * mt.cos(theta[i])) / 100

    # elif param[1] == 125:
    #     if param[0] == 122:
    #         xl = [0, 1.25, 2.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100]
    #         xu = [0, 1.25, 2.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100]
    # elif param[1] == 155:
    # elif param[1] == 185:
    # elif param[1] == 104:
    # elif param[1] == 124:
    # elif param[1] == 144:
    # elif param[1] == 75:
    # elif param[1] == 76:
    # elif param[1] == 126:
    # elif param[1] == 156:
    # elif param[1] == 56:
    # elif param[1] == 176:
    # elif param[1] == 55:
    # elif param[1] == 105:
    # elif param[1] == 175:
    else:
        print('Error: the selected Benzing airfoil does not exist into the database. Please check airfoil number (' + i + ') parameters. \n')
        sys.exit()


    x = (np.concatenate([np.flip(xl, 0), xu[1:]])).tolist()
    y = (np.concatenate([np.flip(yl, 0), yu[1:]])).tolist()

    j = len(y) - 1
    while y[0] > y[j] or abs(y[0] - y[j]) <= 0.001:
        x.pop(0)
        x.pop(len(x) - 1)
        y.pop(0)
        y.pop(len(y) - 1)
        j = len(y) - 1

    if x[0] >= x[len(x) - 1]:
        m = (y[0] - y[1]) / (x[0] - x[1])
        m2 = - 1 / ((y[len(y) - 1] - y[len(y) - 2]) / (x[len(x) - 1] - x[len(x) - 2]))

        q = y[1] - m * x[1]
        q2 = y[len(x) - 1] - m2 * x[len(x) - 1]

        x.pop(0)
        y.pop(0)
        x.insert(0, - (q2 - q) / (m2 - m))
        y.insert(0, m * x[0] + q)
        count = 1
        while count != 0:
            if x[0] <= x[1] and abs(x[0] - x[1]) >= 1e-4:
                x.pop(1)
                y.pop(1)
            elif abs(x[0] - x[1]) < 1e-4:
                x.pop(0)
                y.pop(0)
            else:
                count = 0

    else:
        m = (y[len(y) - 1] - y[len(y) - 2]) / (x[len(x) - 1] - x[len(x) - 2])
        m2 = - 1 / ((y[0] - y[1]) / (x[0] - x[1]))

        q = y[len(y) - 2] - m * x[len(x) - 2]
        q2 = y[0] - m2 * x[0]

        x.pop(len(x) - 1)
        y.pop(len(y) - 1)
        x.append(- (q2 - q) / (m2 - m))
        y.append(m * x[len(x) - 1] + q)
        count = 1
        while count != 0:
            if x[len(x) - 1] <= x[len(x) - 2] and abs(x[len(x) - 1] - x[len(x) - 2]) >= 1e-4:
                x.pop(len(x) - 2)
                y.pop(len(x) - 2)
            elif abs(x[len(x) - 1] - x[len(x) - 2]) < 1e-4:
                x.pop(len(x) - 1)
                y.pop(len(x) - 1)
            else:
                count = 0

    x = np.asarray(x)
    y = np.asarray(y)

    TE = 'Y'

    return x, y, TE
