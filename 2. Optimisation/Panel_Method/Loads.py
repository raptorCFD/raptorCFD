def Loads(p, Cp, alpha):
    import math as mt
    import numpy as np

    alpha = float(alpha)
    n = len(p.panel)

    CFx = 0
    CFy = 0

    for i in range(0, n-1):
        dx = p.panel[i].P2.x - p.panel[i].P1.x
        dy = p.panel[i].P2.y - p.panel[i].P1.y
        CFx = CFx + Cp[i] * dy
        CFy = CFy - Cp[i] * dx

    Cd = CFx * mt.cos(np.deg2rad(alpha)) + CFy * mt.sin(np.deg2rad(alpha))
    Cl = CFx * (-mt.sin(np.deg2rad(alpha))) + CFy * mt.cos(np.deg2rad(alpha))

    return Cl, Cd
