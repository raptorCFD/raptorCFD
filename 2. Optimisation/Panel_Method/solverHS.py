from Panel_Method.AICMatrix import AICMatrix
from Panel_Method.RHSVector import RHSVector
from Panel_Method.Loads import Loads
from Panel_Method.Velocity import Velocity
from multiGeometry import multiGeometry
from Panel_Method.PanelsMulti import PanelsMulti
from Panel_Method.AICMatrixMulti import AICMatrixMulti
from Panel_Method.RHSVectorMulti import RHSVectorMulti
from Panel_Method.VelocityMulti import VelocityMulti
from Panel_Method.Panels import Panels
from Panel_Method.PressureCoeff import PressureCoeff
import numpy as np


def solverHS(alpha, npoint, aname, ref_airfoil, typ, TE, dist, crel, n):

    # Usage:
    # - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha)
    # - [Cl, Cd, Cp, maxdCp] = solverHS(npoint, aname, alpha, dist, crel)

    # maxdCp is a matrix; each row corresponds to an airfoil, first column corresponds to lower part, second column corresponds to upper part
    metaPan = 0

    if len(alpha) == 1:
        p1 = 0
        nairfoils = 1
        xmax = 1
        ymax = 0

        aname1 = aname[0]
        alpha1 = alpha

        # Airfoil discretization and plotting
        x, y, index_Above1, index_Below1, index_Above2, index_Below2, TE_len, ref_airfoil = multiGeometry(nairfoils, npoint, aname, alpha, dist, crel, typ, TE, ref_airfoil)
        p = Panels(x, y)

        # Aerodynamic Influence Coefficients Matrix[AIC]
        AIC = AICMatrix(p)

        # Right Hand Side Vector {RHS}
        RHS = RHSVector(p, alpha1, 1)

        # System solution
        SOL = AIC / RHS

        # Velocity
        v = Velocity(p, alpha1, 1, SOL)  # need this for Cp

        # Pressure
        Cp = 1 - (v ** 2) / 1
        Cpmin = min(Cp)
        maxdCp = abs(Cpmin - Cp[len(Cp) - 1])

        # Aerodynamic coefficients
        Cl, Cd = Loads(p, Cp, alpha1)  # omitted arguments: Cm, CmLE

        # Stuff to achieve uniformity with multiple airfoil case
        # x = {x}
        # y = {y}
        # Cp = {Cp}
        # v = {v}

    elif len(alpha) >= 2:
        nairfoils = len(typ)
        p1 = [0] * nairfoils

        # Load inputs
        # dist = varargin{1}
        # crel = varargin{2}

        # get multi-geometry
        # if count(py.sys.path, '') == 0:
        #    insert(py.sys.path, int32(0), '')

        x, y, index_Above1, index_Below1, index_Above2, index_Below2, TE_len, ref_airfoil = multiGeometry(nairfoils, npoint, aname, alpha, dist, crel, typ, TE, ref_airfoil)

        index = [len(x[ii]) for ii in range(0, nairfoils)]

        # panels
        for i in range(0, nairfoils):  # this runs backwards to avoid preallocation issues!
            np.asarray(x[i])
            np.asarray(y[i])
            p1[i] = Panels(x[i], y[i], index[i])

        p, metaPan = PanelsMulti(p1)

        # Influence matrix
        AIC = AICMatrixMulti(p, metaPan, nairfoils)

        # Right Hand Side Vector {RHS}
        RHS = RHSVectorMulti(p, alpha, 1)

        # Solution
        SOL = np.linalg.solve(AIC, np.asarray(RHS).T)

        # Velocity on profile
        v = VelocityMulti(p, metaPan, nairfoils, alpha, 1, SOL)

        #print('v1', v[1])

        # Preallocation of coefficients
        Cp, maxdCp, Cl, Cd, xmax, ymax = [0] * nairfoils, [0] * nairfoils, [0] * nairfoils, [0] * nairfoils, [0] * nairfoils, [0] * nairfoils

        # Calculation of coefficients
        for i in range(0, nairfoils):
            xmax[i] = max(x[i])
            ymax[i] = max(y[i])
            Cp[i] = PressureCoeff(v[i], 1)
            Cl[i], Cd[i] = Loads(p1[i], Cp[i], alpha[i])
            # calculation of maxdCp for Valarezo - Chin
            minCp = min(Cp[i])
            maxdCp[i] = abs(minCp - Cp[i][len(Cp[i]) - 1])

    else:
        import sys
        print('wrong input; see documentation for instructions on how to use this function.')
        sys.exit()

    return Cl, Cd, xmax, ymax, maxdCp, x, y, index, p, p1, SOL, metaPan, nairfoils
