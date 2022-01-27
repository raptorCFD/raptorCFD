from Panel_Method.ConstantSource2D import ConstantSource2D
from Panel_Method.Rotation import Rotation
import numpy as np


def AICMatrixMulti(p, metaPan, nairfoils):
    import math as mt

    # extract input
    npan = metaPan.npan
    idx_zeroPan = metaPan.idx_zeroPan

    # preallocation
    ntot = len(p.panel)

    AIC = [None] * ntot

    for i in range(0, ntot):
        AIC_new = [0] * ntot

        ni = np.asarray([-mt.sin(p.panel[i].beta), mt.cos(p.panel[i].beta)])

        # sotto matrice Aij(sorgenti su sorgenti) e bs(sorgenti su vortici)
        for j in range(0, ntot - nairfoils):

            if i == j:
                Rt = np.asarray(Rotation(p.panel[j].beta)).T
                us = np.asarray([[0], [0.5]])
                us = Rt.dot(us)

                AIC_new[j] = ni.dot(us)[0]

            else:
                us = np.asarray(ConstantSource2D(1, p.panel[j], p.panel[i]))
                us = np.asarray([[us[0]], [us[1]]])

            AIC_new[j] = ni.dot(us)[0]

       # sottomatrice b_v(vortici su sorgenti)
        for k in range(0, nairfoils):

            for j in range(idx_zeroPan[k], npan[k] + idx_zeroPan[k]):

                if i == j:
                    Rt = np.asarray(Rotation(p.panel[j].beta)).T
                    uv = Rt.dot(np.asarray([[0.5], [0]]))

                else:

                    us = ConstantSource2D(1, p.panel[j], p.panel[i])
                    uv = np.asarray([[us[1]], [-us[0]]])

                AIC_new[ntot - nairfoils + k] = AIC_new[ntot - nairfoils + k] + ni.dot(uv)[0]

        AIC[i] = AIC_new

    return AIC
