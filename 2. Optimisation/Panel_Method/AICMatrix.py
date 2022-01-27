from Panel_Method.Rotation import Rotation
from Panel_Method.ConstantSource2D import ConstantSource2D

def AICMatrix(p):
    import math as mt
    import numpy as np
    n = len(p.panel)
    
    AIC = [[0] * n] * n
    
    for i in range(0, n):
    
        ni = np.asarray([-mt.sin(p.panel(i).beta), mt.cos(p.panel(i).beta)])
        
        for j in range(1, n-1):
        
            if i == j:
                Rt = np.asarray(Rotation(p.panel(j).beta)).T
                us = np.asarray([[0], [0.5]])
                uv = np.asarray([[0.5], [0]])
                us = Rt.dot(us)
                uv = Rt.dot(uv)
            else:
                us = ConstantSource2D(1, p.panel(j), p.panel(i))
                uv = np.asarray([[us[1]], [-us[0]]])

            AIC[i][j] = ni.dot(us)
            AIC[i][n] = AIC[i][n] + ni.dot(uv)
    
    return AIC
