from Panel_Method.ConstantSource2D import ConstantSource2D_v
from Panel_Method.ConstantVortex2D import ConstantVortex2D_v


def Velocity(p, alpha, U, SOL):
    import math as mt
    import numpy as np
    n = len(p.panel)
    uInf = U*[mt.cos(np.deg2rad(alpha)), mt.sin(np.deg2rad(alpha))]
    
    v = [0] * (n-1)

    ti = [mt.cos([p.panel[:(n-1)].beta]), mt.sin([p.panel[:(n-1)].beta])]
    
    # preallocation
    us, uv = [[0] * (n-1)] * 2
    idx_old = np.asarray(range(0, n-1)).tolist()

    for j in range(0, n-1):
      idx = idx_old
      idx.pop(j)

      us_new = ConstantSource2D_v(SOL[j], p.panel[j], p.panel[idx])
      uv_new = ConstantVortex2D_v(SOL[len(SOL)-1], p.panel[j], p.panel[idx])

      for k in range(0, n-1):
          if k != j:
              us[:, k] = us_new[:, k]
              uv[:, k] = uv_new[:, k]
          else:
              us[:, j] = p.panel[j].R * [0, 0.5] * SOL[j]
              uv[:, j] = p.panel[j].R * [0.5, 0] * SOL[len(SOL) - 1]

      us[:,idx] = ConstantSource2D_v(SOL[j], p.panel[j], p.panel(idx))
      uv[:,idx] = ConstantVortex2D_v(SOL[len(SOL)-1], p.panel[j], p.panel(idx))


      v = v + sum(ti*us,1) + sum(ti*uv,1)

    v = v + ti*uInf
    
    return v