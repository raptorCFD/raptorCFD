def RHSVectorMulti(p, alpha, U):
    import math as mt
    import numpy as np
    n = len(p.panel)

    alpha1 = float(alpha[0]) # il sistema di riferimento globale è solidale al primo profilo!

    uInf1 = [U * mt.cos(np.deg2rad(alpha1)), U * mt.sin(np.deg2rad(alpha1))]

    RHS = [0] * n

    for i in range(0, n):
        RHS[i] = mt.sin(p.panel[i].beta) * uInf1[0] - mt.cos(p.panel[i].beta) * uInf1[1]

    return RHS