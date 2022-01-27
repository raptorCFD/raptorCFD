import math as mt


def wall_spacing(Re, crel, mu, nairfoils, y_plus, progr, thick, temp, TE, Mach):
    U = Mach * mt.sqrt(1.4 * 287 * temp)
    nu = U / Re
    rho = mu / nu
    Rex, cf, tau_w, u_norm, s, BL_thickness, norm_nodes, norm_nodes_TE = [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [
        None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [None] * nairfoils, [None] * nairfoils

    for j in range(0, nairfoils):
        Rex[j] = Re * crel[j]
        # Boundary Layer thickness
        # https://en.wikipedia.org/wiki/Boundary_layer_thickness#Momentum_Thickness
        # BL_thickness[j] = 0.37 * crel[j] / (Rex[j] ** 0.2)
        BL_thickness[j] = 0.37 * crel[j] / (Rex[j] ** 0.2) * thick[j]
        cf[j] = 0.0576 * Rex[j] ** (-0.2)
        tau_w[j] = 0.5 * cf[j] * rho * U ** 2
        u_norm[j] = mt.sqrt(tau_w[j] / rho)
        s[j] = y_plus[j] * nu / (rho * u_norm[j])
        norm_nodes[j] = int(mt.log(1 + (progr[j] - 1) * BL_thickness[j] / s[j], progr[j]))  # resolution of geometric sequence
        if TE[j] == 'Y':
            norm_nodes_TE[j] = int(mt.log(1 + (progr[j] - 1) * (BL_thickness[j] * 0.01) / s[j], progr[j]))
        else:
            pass

    return BL_thickness, norm_nodes, s, progr, Rex, norm_nodes_TE, U, rho
