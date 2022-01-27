def velocityField(p, nairfoils, alpha, U, SOL, xa, ya, metaPan):
    import math as mt
    import numpy as np

    def ConstantSource2DField(sigma_J, panel_J, x, y):

        x1 = -panel_J.d / 2
        x2 = panel_J.d / 2

        x_transpose = [x[i][0] for i in range(0, len(x))]
        y_transpose = [y[i][0] for i in range(0, len(y))]

        P_I = (np.asarray(panel_J.R).dot([[x_transpose[i] - panel_J.C.x for i in range(0, len(x_transpose))], [y_transpose[i] - panel_J.C.y for i in range(0, len(y_transpose))]])).tolist()

        upx, upy = [0] * len(x), [0] * len(x)
        # Local
        for i in range(0, len(upx)):
            upx[i] = (sigma_J / (4 * mt.pi)) * mt.log(((P_I[0][i] - x1) ** 2 + P_I[1][i] ** 2) / ((P_I[0][i] - x2) ** 2 + P_I[1][i] ** 2))
            upy[i] = (sigma_J / (2 * mt.pi)) * (mt.atan((P_I[0][i] - x1) / P_I[1][i]) - mt.atan((P_I[0][i] - x2) / P_I[1][i]))

        UU = np.asarray(panel_J.R).T.dot(np.asarray([upx, upy]))
        u = UU[0, :]
        v = UU[1, :]

        return u, v

    def ConstantVortex2DField(gamma_J, panel_J, x, y):
        # R = Rotation (panel_J.beta);

        x1 = -panel_J.d / 2
        x2 = panel_J.d / 2

        x_transpose = [x[i][0] for i in range(0, len(x))]
        y_transpose = [y[i][0] for i in range(0, len(y))]

        P_I = (np.asarray(panel_J.R).dot([[x_transpose[i] - panel_J.C.x for i in range(0, len(x_transpose))],
                                          [y_transpose[i] - panel_J.C.y for i in range(0, len(y_transpose))]])).tolist()

        upx, upy = [0] * len(x), [0] * len(x)

        for i in range(0, len(upx)):
            upx[i] = (gamma_J / (2 * mt.pi)) * (mt.atan((P_I[0][i] - x1) / P_I[1][i]) - mt.atan((P_I[0][i] - x2) / P_I[1][i]))
            upy[i] = -(gamma_J / (4 * mt.pi)) * mt.log(((P_I[0][i] - x1) ** 2 + P_I[1][i] ** 2) / ((P_I[0][i] - x2) ** 2 + P_I[1][i] ** 2))

        UU = np.asarray(panel_J.R).T.dot(np.asarray([upx, upy]))
        u = UU[0, :]
        v = UU[1, :]

        return u, v

    if nairfoils == 1:
        npan = len(p.panel) - 1
        idx_zeroPan = 0
    else:
        npan = metaPan.npan
        idx_zeroPan = metaPan.idx_zeroPan


    ntot = len(p.panel)

    alpha1 = float(alpha[0])  # s.o.f. referred to the first element of the configuration
    uInf1 = [U * mt.cos(np.deg2rad(alpha1)), U * mt.sin(np.deg2rad(alpha1))]


    # get boundaries for domain
    xmin = min([p.panel[i].C.x for i in range(0, ntot)]) - 0.75
    xmax = max([p.panel[i].C.x for i in range(0, ntot)]) + 0.75
    ymin = min([p.panel[i].C.y for i in range(0, ntot)]) - 0.75
    ymax = max([p.panel[i].C.y for i in range(0, ntot)]) + 0.75

    # build domain
    Ngrid = 150
    xgrid = np.linspace(xmin, xmax, Ngrid).T
    ygrid = np.linspace(ymin, ymax, Ngrid).T

    [Xgrid, Ygrid] = np.meshgrid(xgrid, ygrid)

    xpoints, ypoints = [0] * (Ngrid ** 2), [0] * (Ngrid ** 2)

    s = 0
    for j in range(0, Ngrid):
        for i in range(0, Ngrid):
            xpoints[s] = [Xgrid[i][j]]
            ypoints[s] = [Ygrid[i][j]]
            s += 1

    # preallocate velocity components
    u = np.asarray([0] * (Ngrid ** 2)).T
    v = np.asarray([0] * (Ngrid ** 2)).T

    for j in range(0, ntot - nairfoils):
        us, vs = ConstantSource2DField(SOL[j], p.panel[j], xpoints, ypoints)
        u = u + us
        v = v + vs

    for k in range(0, nairfoils):
        for j in range(0, npan[k] + idx_zeroPan[k]):  # ((k - 1) * nterz + 1): k * nterz
            uv, vv = ConstantVortex2DField(SOL[len(SOL) - 1 - nairfoils + k], p.panel[j], xpoints, ypoints)
            u = u + uv
            v = v + vv

    u = u + uInf1[0]
    v = v + uInf1[1]

    from matplotlib import path
    for k in range(0, nairfoils):
        coord = [None] * len(xa[k])
        for kk in range(0, len(coord)):
            coord[kk] = [xa[k], ya[k]]
        p = path.Path(coord)

        for ii in range(0, len(xpoints)):
            if p.contains_points([(xpoints[ii], ypoints[ii])]):
                u[ii] = 0
                v[ii] = 0

    u = np.reshape(u, (Ngrid, Ngrid))
    v = np.reshape(v, (Ngrid, Ngrid))

    return u, v, Xgrid, Ygrid