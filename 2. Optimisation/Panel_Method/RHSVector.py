def RHSVector (p, alpha, U):
    import math as mt
    import numpy as np
    n = len(p.panel)
    uInf = [U*mt.cos(np.deg2rad(alpha)), U*mt.sin(np.deg2rad(alpha))]
    RHS = [0] * n

    for i in range(0, n):
        RHS[i] = mt.sin(p.panel[i].beta) * uInf[0] - mt.cos(p.panel[i].beta) * uInf[1]

    return RHS
