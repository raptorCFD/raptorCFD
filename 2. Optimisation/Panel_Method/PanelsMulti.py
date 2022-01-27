class ClassMetapan:
    def __init__(self, npan, idx_zeroPan):
        self.npan = npan
        self.idx_zeroPan = idx_zeroPan


class PanelClass:
    def __int__(self, panel):
        self.panel = panel


class PanelDetails:
    def __init__(self, P1, P2, C, d, beta, R):
        self.P1 = P1
        self.P2 = P2
        self.C = C
        self.d = d
        self.beta = beta
        self.R = R


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def PanelsMulti(p1):
    nairfoil = len(p1)

    # get number of pts for each airfoil
    nptz, npan = [0] * nairfoil, [0] * nairfoil
    for ii in range(0, nairfoil):
        nptz[ii] = len(p1[ii].panel)
        # get number of TRUE panels(excludes Kutta panels) for each airfoil
        npan[ii] = nptz[ii] - 1

    # zero index for panels of each airfoil
    idx_zeroPan = [0] * nairfoil

    for ii in range(1, nairfoil):
        idx_zeroPan[ii] = sum(npan[:ii])

    # zero index for kutta panels
    idx_zeroKutta = sum(npan)

    p = PanelClass()
    p.panel = [0] * (sum(npan) + nairfoil)

    # ntot = sum(nptzsa)

    # nterz = (ntot - nairfoil) / nairfoil

    for k in range(0, nairfoil):

        # nterz = npoints - 1

        for j in range(0, npan[k]):

            i = idx_zeroPan[k] + j

            P1 = Point(p1[k].panel[j].P1.x, p1[k].panel[j].P1.y)
            P2 = Point(p1[k].panel[j].P2.x, p1[k].panel[j].P2.y)
            C = Point(p1[k].panel[j].C.x, p1[k].panel[j].C.y)
            d = p1[k].panel[j].d
            beta = p1[k].panel[j].beta
            R = p1[k].panel[j].R
            p.panel[i] = PanelDetails(P1, P2, C, d, beta, R)


        i = idx_zeroKutta + k

        npoints = nptz[k]

        P1 = Point(p1[k].panel[npoints-1].P1.x, p1[k].panel[npoints-1].P1.y)
        P2 = Point(p1[k].panel[npoints-1].P2.x, p1[k].panel[npoints-1].P2.y)
        C = Point(p1[k].panel[npoints-1].C.x, p1[k].panel[npoints-1].C.y)

        d = p1[k].panel[npoints-1].d
        beta = p1[k].panel[npoints-1].beta
        R = p1[k].panel[npoints-1].R
        p.panel[i] = PanelDetails(P1, P2, C, d, beta, R)

    metaPan = ClassMetapan
    metaPan.npan = npan
    metaPan.idx_zeroPan = idx_zeroPan

    return p, metaPan
