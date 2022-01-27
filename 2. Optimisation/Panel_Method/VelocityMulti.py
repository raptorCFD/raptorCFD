from Panel_Method.Rotation import Rotation
from Panel_Method.ConstantSource2D import ConstantSource2D
from Panel_Method.ConstantVortex2D import ConstantVortex2D


def VelocityMulti(p, metaPan, nairfoils, alpha, U, SOL):
    import numpy as np
    import math as mt

    # extract input
    npan = metaPan.npan
    idx_zeroPan = metaPan.idx_zeroPan
    
    ntot = len(p.panel)
    
    alpha1 = float(alpha[0])    # il sist.di riferimento è solidale al primo profilo
    uInf1 = [U * mt.cos(np.deg2rad(alpha1)), U * mt.sin(np.deg2rad(alpha1))]
    
    v = [0] * (ntot - nairfoils)

    for i in range(0, ntot - nairfoils):
    
        ti = np.asarray([mt.cos(p.panel[i].beta), mt.sin(p.panel[i].beta)])

        # CONTRIBUTO SORGENTI

        for j in range(0, ntot - nairfoils):

            if i == j:
                Rt = np.asarray(Rotation(p.panel[j].beta)).T
                us = (Rt.dot([[0], [0.5]])).dot(SOL[j])
                # uv = Rt * [0.5, 0]*SOL(end)
            else:
                us = ConstantSource2D(SOL[j], p.panel[j], p.panel[i])

            v[i] = v[i] + ti.dot(us) #+ ti * uv

        # CONTRIBUTO VORTICI

        for k in range(0, nairfoils):

            for j in range(idx_zeroPan[k], npan[k] + idx_zeroPan[k]):

                if i == j:
                    Rt = np.asarray(Rotation(p.panel[j].beta)).T
                    uv = (Rt.dot([[0.5], [0]]))*(SOL[len(SOL) - nairfoils + k])
                else:
                    uv = ConstantVortex2D(SOL[len(SOL) - nairfoils + k], p.panel[j], p.panel[i])

                v[i] = v[i] + ti.dot(uv)

        v[i] = (v[i] + ti.dot(uInf1)).tolist()[0]

    vaux = [0] * nairfoils

    for k in range(0, nairfoils):
        vaux[k] = [0] * npan[k]
        for i in range(0, npan[k]):
            if k == 0:
                vaux[k][i] = v[i]
            else:
                vaux[k][i] = v[i + idx_zeroPan[k]]

    return vaux


