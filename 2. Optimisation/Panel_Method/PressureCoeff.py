def PressureCoeff(v, U):

    Cp = [0] * len(v)

    for i in range(0, len(v)):
        Cp[i] = 1-(v[i]**2)/(U**2)

    return Cp
