import math as mt
from Panel_Method.Rotation import Rotation


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


def Panels(x, y, n=0):
    p = PanelClass()
    p.panel = [0] * n

    for i in range(0, n-1):
        
        d = mt.sqrt((x[i+1] - x[i]) ** 2 + (y[i+1] - y[i]) ** 2)
        beta = mt.atan2((y[i+1] - y[i]), (x[i+1] - x[i]))
        R = Rotation(beta)
        P1 = Point(x[i], y[i])
        P2 = Point(x[i+1], y[i+1])
        C = Point((x[i] + x[i+1]) / 2, (y[i] + y[i+1]) / 2)
        p.panel[i] = PanelDetails(P1, P2, C, d, beta, R)

    # Add Wake Panel
    panel_1 = p.panel[0]
    panel_2 = p.panel[n - 2]
    P1 = Point(panel_1.P1.x, panel_1.P1.y)
    P2 = Point((panel_1.P1.x + (panel_1.P1.x - panel_1.P2.x) + panel_2.P2.x + (panel_2.P2.x - panel_2.P1.x)) / 2, (panel_1.P1.y + (panel_1.P1.y - panel_1.P2.y) + panel_2.P2.y + (panel_2.P2.y - panel_2.P1.y)) / 2)

    d = mt.sqrt((P2.x - P1.x) ** 2 + (P2.y - P1.y) ** 2)
    C = Point((P1.x + P2.x) / 2, (P1.y + P2.y) / 2)
    beta = mt.atan2(P2.y - P1.y, P2.x - P1.x)

    # Rotation
    R = Rotation(beta)
    p.panel[n-1] = PanelDetails(P1, P2, C, d, beta, R)

    return p
