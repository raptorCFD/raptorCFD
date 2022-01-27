#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math as mt
import sys


def camberlineIGP(c1, c2, c3, c4, k):
    xc = (3 * c1 * k * (1 - k) ** 2) + (3 * c2 * (1 - k) * k ** 2) + k ** 3
    yc = (3 * c3 * k * (1 - k) ** 2) + (3 * c4 * (1 - k) * k ** 2)

    return xc, yc


def camberlineNACA4(xc, m, p, t, c):
    yt = 5 * t * (0.2969 * mt.sqrt(xc) - 0.126 * xc - 0.3516 * (xc ** 2) + 0.2843 * (xc ** 3) - 0.1015 * (
            xc ** 4))

    if p == 0:
        yc = 0
        theta = mt.atan(0)
    elif xc <= p * c:
        yc = m / (p ** 2) * (2 * p * (xc / c) - (xc / c) ** 2)
        theta = mt.atan(2 * m / (p ** 2) * (p - xc / c))
    else:
        yc = m / ((1 - p) ** 2) * ((1 - 2 * p) + 2 * p * (xc / c) - (xc / c) ** 2)
        theta = mt.atan(2 * m / ((1 - p) ** 2) * (p - xc / c))

    return yc, theta, yt


def camberlineNACA5_NONREFL(xc, c, param):
    yt = 5 * param[3]/100 * (0.2969 * mt.sqrt(xc) - 0.126 * xc - 0.3516 * (xc ** 2) + 0.2843 * (xc ** 3) - 0.1015 * (
            xc ** 4))
    if int(str(param[0]) + str(param[1]) + str(param[2])) == 210:
        p = 0.05
        m = 0.0580
        k1 = 361.40
    elif int(str(param[0]) + str(param[1]) + str(param[2])) == 220:
        p = 0.10
        m = 0.126
        k1 = 51.640
    elif int(str(param[0]) + str(param[1]) + str(param[2])) == 230:
        p = 0.15
        m = 0.2025
        k1 = 15.957
    elif int(str(param[0]) + str(param[1]) + str(param[2])) == 240:
        p = 0.20
        m = 0.290
        k1 = 6.643
    elif int(str(param[0]) + str(param[1]) + str(param[2])) == 250:
        p = 0.25
        m = 0.391
        k1 = 3.230
    else:
        print("Error:" + param[0] + param[1] + param[2] + " is not included in database")
        sys.exit()

    if xc <= p * c:
        yc = k1 / 6 * (((xc / c) ** 3) - 3 * m * ((xc / c) ** 2) + (m ** 2) * (3 - m) * xc / c)
        theta = (k1/6) * ((3 * ((xc / c) ** 2)) - (6 * m * (xc / c)) + (m ** 2) * (3 - m))
    else:
        yc = ((k1 * (m ** 3)) / 6) * (1 - (xc / c))
        theta = (-1 * k1 * (m ** 3)) / 6

    return yc, theta, yt


def camberlineNACA5_REFL(xc, c, param):
    yt = 5 * param[3]/100 * (0.2969 * mt.sqrt(xc) - 0.126 * xc - 0.3516 * (xc ** 2) + 0.2843 * (xc ** 3) - 0.1015 * (
            xc ** 4))

    if int(str(param[0]) + str(param[1]) + str(param[2])) == 221:
        p = 0.10
        m = 0.130
        k1 = 51.990
        k2_k1 = 0.000764
    elif int(str(param[0]) + str(param[1]) + str(param[2])) == 231:
        p = 0.15
        m = 0.217
        k1 = 15.793
        k2_k1 = 0.00677
    elif int(str(param[0]) + str(param[1]) + str(param[2])) == 241:
        p = 0.20
        m = 0.318
        k1 = 6.520
        k2_k1 = 0.0303
    elif int(str(param[0]) + str(param[1]) + str(param[2])) == 251:
        p = 0.25
        m = 0.441
        k1 = 3.191
        k2_k1 = 0.1355
    else:
        print("Error:" + param[0] + param[1] + param[2] + " is not included in database")
        sys.exit()

    if xc <= m * c:
        yc = k1 / 6 * (((xc/c - m)**3) - k2_k1 * ((1-m)**3) * xc /c - (m**3) * xc / c + (m**3))
        theta = k1 / 6 * (3 * ((xc/c - m) ** 2) - k2_k1 * ((1 - m)**3) - (m**3))
    else:
        yc = k1 / 6 * (k2_k1 * ((xc / c - m) ** 3) - k2_k1 * ((1 - m) ** 3) * xc / c - (m ** 3) * xc / c + (m ** 3))
        theta = k1 / 6 * (3 * k2_k1 * ((xc/c - m) ** 2) - k2_k1 * ((1 - m)**3) - (m**3))
    return yc, theta, yt


def camberlineNACA4_MOD(xc, c, m, p, t, I, T):

    if I <= 8:
        r_le = c * 1.1019 * ((I * t / (6 * c)) ** 2)
        Xi = I / 6
    elif I == 9:
        r_le = 3 * c * 1.1019 * ((t / c) ** 2)
        Xi = 10.3933
    else:
        print("Error: fourth parameter of NACA 4-digit modified not found. Check the value inserted. \n\n")
        sys.exit()

    # Camberline is the same of NACA 4-digit
    if p == 0:
        yc = 0
        theta = mt.atan(0)
    elif xc <= p * c:
        yc = m / (p ** 2) * (2 * p * (xc / c) - (xc / c) ** 2)
        theta = mt.atan(2 * m / (p ** 2) * (p - xc / c))
    else:
        yc = m / ((1 - p) ** 2) * ((1 - 2 * p) + 2 * p * (xc / c) - (xc / c) ** 2)
        theta = mt.atan(2 * m / ((1 - p) ** 2) * (p - xc / c))

    # d1 determination
    if T == 0.2:
        d1 = 0.200
    elif T == 0.3:
        d1 = 0.234
    elif T == 0.4:
        d1 = 0.315
    elif T == 0.5:
        d1 = 0.465
    elif T == 0.6:
        d1 = 0.700
    else:
        print('Error: last value inserted about NACA 4-digit modified not found.')
        sys.exit()

    d2 = (0.294 - 2 * (1-T)*d1) / ((1-T)**2)
    d3 = (-0.196 + (1-T)*d1) / ((1-T)**3)

    a0 = 0.296904 * Xi

    rho1 = (1/5) * (((1-T)**2)/(0.588 - 2*d1*(1-T)))

    a1 = (0.3 / T) - (15 * a0 / (8 * mt.sqrt(T))) - (T / (10 * rho1))

    a2 = -(0.3 / (T ** 2)) + (5 * a0 / (4 * (T ** (3 / 2)))) + (1 / (5 * rho1))

    a3 = (0.1 / (T**3)) - (0.375 * a0 / (T ** (5/2))) - (1 / (10 * rho1 * T))

    if xc <= T * c:
        yt = 5 * t * (a0 * mt.sqrt(xc/c) + a1 * (xc/c) + a2 * (xc/c) ** 2 + a3 * (xc/c) ** 3)
    else:
        yt = 5 * t * (0.002 + d1 * (1 - xc/c) + d2 * (1 - xc/c) ** 2 + d3 * (1 - xc/c) ** 3)

    return yc, theta, yt



