import math as mt
import numpy as np


def ConstantVortex2D(gamma_J, panel_J, panel_I):
    #from Rotation import Rotation
    #
    #R = Rotation (panel_J.beta)
    
    x1 = -panel_J.d/2
    x2 = panel_J.d/2
    
    P_I = (np.asarray(panel_J.R).dot(np.asarray([panel_I.C.x-panel_J.C.x, panel_I.C.y-panel_J.C.y]))).tolist()
    upx = (gamma_J / (2 * mt.pi)) * (mt.atan((P_I[0] - x1) / P_I[1]) - mt.atan((P_I[0] - x2) / P_I[1]))
    upy = -(gamma_J / (4 * mt.pi)) * mt.log(((P_I[0] - x1) **2 + P_I[1] **2) / ((P_I[0] - x2) **2 + P_I[1] **2))
    
    u = (np.asarray(panel_J.R).T.dot(np.asarray([upx, upy]))).tolist()

    return u


def ConstantVortex2D_v(gamma_J, panel_J, panels):
    x1 = -panel_J.d/2
    x2 = panel_J.d/2

    CPan = [panels.C]

    P_I = (np.asarray(panel_J.R).dot(np.asarray([[CPan.x]-panel_J.C.x, [CPan.y]-panel_J.C.y]))).tolist()

    upx = (gamma_J/(2*mt.pi)) * (mt.atan((P_I[0]-x1)/P_I[1]) - mt.atan((P_I[0]-x2)/P_I[1]))
    upy = -(gamma_J/(4*mt.pi)) * mt.log(((P_I[0]-x1)**2 + P_I[1]**2)/((P_I[0]-x2)**2 + P_I[1]**2))

    u = (np.asarray(panel_J.R).dot(np.asarray([upx, upy]))).tolist()

    return u
