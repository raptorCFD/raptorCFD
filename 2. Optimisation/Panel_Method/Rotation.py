def Rotation(beta):
    import math as mt
    R = [[mt.cos(beta), mt.sin(beta)], [-mt.sin(beta), mt.cos(beta)]]
    return R
