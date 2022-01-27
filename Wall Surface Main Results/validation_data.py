import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# function to convert to superscript
def get_super(x):
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    super_s = "ᴬᴮᶜᴰᴱᶠᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾQᴿˢᵀᵁⱽᵂˣʸᶻᵃᵇᶜᵈᵉᶠᵍʰᶦʲᵏˡᵐⁿᵒᵖ۹ʳˢᵗᵘᵛʷˣʸᶻ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾"
    res = x.maketrans(''.join(normal), ''.join(super_s))
    return x.translate(res)

# function to convert to subscript
def get_sub(x):
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    sub_s = "ₐ₈CDₑբGₕᵢⱼₖₗₘₙₒₚQᵣₛₜᵤᵥwₓᵧZₐ♭꜀ᑯₑբ₉ₕᵢⱼₖₗₘₙₒₚ૧ᵣₛₜᵤᵥwₓᵧ₂₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎"
    res = x.maketrans(''.join(normal), ''.join(sub_s))
    return x.translate(res)


def validation_data(MY_FILE, EXP_DATA_CP, EXP_DATA_CF, my_label, exp_label, my_size, exp_size, colour_exp_cp, colour_exp_cf):
    df = pd.read_csv(MY_FILE)

    my_x = df["x"].tolist()
    my_y = df["y"].tolist()

    try:
        my_rho = df["Density"].tolist()
        print('The maximum value of wall density is (' + str(round(max(my_rho), 3)) + ') at ' + str(
            round(my_x[np.argmax(my_rho)], 3)) + ' m.\n')
        print('The minimum value of wall density is (' + str(round(min(my_rho), 3)) + ') at ' + str(
            round(my_x[np.argmin(my_rho)], 3)) + ' m.\n\n')
        plt.scatter(my_x, my_rho, c='black', s=my_size, label=my_label)
        plt.suptitle('Density(x)')
        plt.ylabel('Density [$kg/m^3$]')
        plt.xlabel('x [m]')
        plt.grid(True)
        if min(my_rho) >= 0:
            min_my_rho = min(my_rho) * 0.9
        else:
            min_my_rho = min(my_rho) * 1.1
        if max(my_rho) >= 0:
            max_my_rho = max(my_rho) * 1.1
        else:
            max_my_rho = max(my_rho) * 0.9

        plt.ylim((min_my_rho, max_my_rho))
        plt.xlim((min(my_x) - 0.1, max(my_x) * 1.1))
        plt.legend()
        plt.savefig('Density(x).png')
        plt.show()
        plt.close()

    except:
        print("Density not available in provided CSV file.\nRespective plot will be skipped.")

    try:
        my_energy = df["Energy"].tolist()
        print('The maximum value of enrrgy is (' + str(round(max(my_energy), 3)) + ') at ' + str(
            round(my_x[np.argmax(my_energy)], 3)) + ' m.\n')
        print('The minimum value of energy is (' + str(round(min(my_energy), 3)) + ') at ' + str(
            round(my_x[np.argmin(my_energy)], 3)) + ' m.\n\n')
        plt.scatter(my_x, my_energy, c='black', s=my_size, label=my_label)
        plt.suptitle('Energy(x)')
        plt.ylabel('Energy')
        plt.xlabel('x [m]')
        plt.grid(True)
        if min(my_energy) >= 0:
            min_my_energy = min(my_energy) * 0.95
        else:
            min_my_energy = min(my_energy) * 1.05
        if max(my_energy) >= 0:
            max_my_energy = max(my_energy) * 1.05
        else:
            max_my_energy = max(my_energy) * 0.95

        plt.ylim((min_my_energy, max_my_energy))
        plt.xlim((min(my_x) - 0.1, max(my_x) * 1.1))
        plt.legend()
        plt.savefig('Energy(x).png')
        plt.show()
        plt.close()
    except:
        print("Energy not available in provided CSV file.\nRespective plot will be skipped.")

    try:
        my_pressure = df["Pressure"].tolist()
        print('The maximum value of wall pressure is (' + str(round(max(my_pressure), 3)) + ') at ' + str(round(my_x[np.argmax(my_pressure)], 3)) + ' m.\n')
        print('The minimum value of wall pressure is (' + str(round(min(my_pressure), 3)) + ') at ' + str(round(my_x[np.argmin(my_pressure)], 3)) + ' m.\n\n')
        plt.scatter(my_x, my_pressure, c='black', s=my_size, label=my_label)
        plt.suptitle('Pressure(x)')
        plt.ylabel('Pressure [K]')
        plt.xlabel('x [m]')
        plt.grid(True)
        if min(my_pressure) >= 0:
            min_my_pressure = min(my_pressure) * 0.95
        else:
            min_my_pressure = min(my_pressure) * 1.05
        if max(my_pressure) >= 0:
            max_my_pressure = max(my_pressure) * 1.05
        else:
            max_my_pressure = max(my_pressure) * 0.95

        plt.ylim((min_my_pressure, max_my_pressure))
        plt.xlim((min(my_x) - 0.1, max(my_x) * 1.1))
        plt.legend()
        plt.savefig('Pressure(x).png')
        plt.show()
        plt.close()
    except:
        print("Pressure not available in provided CSV file.\nRespective plot will be skipped.")

    try:
        my_temp = df["Temperature"].tolist()
        print('The maximum value of wall temperature is (' + str(round(max(my_temp), 3)) + ') at ' + str(
            round(my_x[np.argmax(my_temp)], 3)) + ' m.\n')
        print('The minimum value of wall temperature is (' + str(round(min(my_temp), 3)) + ') at ' + str(
            round(my_x[np.argmin(my_temp)], 3)) + ' m.\n\n')
        plt.scatter(my_x, my_temp, c='black', s=my_size, label=my_label)
        plt.suptitle('Temperature(x)')
        plt.ylabel('Temperature [K]')
        plt.xlabel('x [m]')
        plt.grid(True)
        if min(my_temp) >= 0:
            min_my_temp = min(my_temp) * 0.98
        else:
            min_my_temp = min(my_temp) * 1.02
        if max(my_temp) >= 0:
            max_my_temp = max(my_temp) * 1.02
        else:
            max_my_temp = max(my_temp) * 0.98

        plt.ylim((min_my_temp, max_my_temp))
        plt.xlim((min(my_x) - 0.1, max(my_x) * 1.1))
        plt.savefig('Temperature(x).png')
        plt.legend()
        plt.show()
        plt.close()
    except:
        print("Temperature not available in provided CSV file.\nRespective plot will be skipped.")

    try:
        my_cp = df["Pressure_Coefficient"].tolist()
        print('The maximum value of wall pressure is (' + str(round(max(my_cp), 5)) + ') at ' + str(
            round(my_x[np.argmax(my_cp)], 3)) + ' m.\n')
        print('The minimum value of wall pressure is (' + str(round(min(my_cp), 5)) + ') at ' + str(
            round(my_x[np.argmin(my_cp)], 3)) + ' m.\n\n')
        plt.scatter(my_x, my_cp, c='black', s=my_size, label=my_label)
        plt.suptitle('Pressure coefficient(x)')
        plt.ylabel('$C_P$ [K]')
        plt.xlabel('x [m]')
        plt.grid(True)
        if min(my_cp) >= 0:
            min_my_cp = min(my_cp) * 0.9
        else:
            min_my_cp = min(my_cp) * 1.1
        if max(my_cp) >= 0:
            max_my_cp = max(my_cp) * 1.1
        else:
            max_my_cp = max(my_cp) * 0.9

        plt.ylim((min_my_cp, max_my_cp))
        plt.xlim((min(my_x) - 0.1, max(my_x) * 1.1))
        plt.legend()
        plt.savefig('Cp(x).png')
        plt.show()
        plt.close()

        if EXP_DATA_CP is not None:
            with open(EXP_DATA_CP) as f:
                exp_data = [[float(row.split()[0]), float(row.split()[1])] for row in f]
                exp_data_x = [exp_data[i][0] for i in range(0, len(exp_data))]
                exp_data_cp = [exp_data[i][1] for i in range(0, len(exp_data))]

                plt.scatter(my_x, my_cp, c='black', s=my_size, label=my_label)
                plt.scatter(exp_data_x, exp_data_cp, c=colour_exp_cp, s=exp_size, label=exp_label)
                plt.suptitle('Pressure Coefficient(x)')
                plt.ylabel('Pressure Coefficient []')
                plt.xlabel('x [m]')
                plt.grid(True)
                if min(min(my_cp), min(exp_data_cp)) >= 0:
                    min_cp = min(min(my_cp), min(exp_data_cp)) * 0.95
                else:
                    min_cp = min(min(my_cp), min(exp_data_cp)) * 1.05
                if max(max(my_cp), max(exp_data_cp)) >= 0:
                    max_cp = max(max(my_cp), max(exp_data_cp)) * 1.05
                else:
                    max_cp = max(max(my_cp), max(exp_data_cp)) * 0.95

                plt.ylim((min_cp, max_cp))
                plt.xlim((min(my_x) - 0.1, max(max(my_x), max(exp_data_x)) * 1.1))
                plt.legend()
                plt.savefig('Cp(x): comparison.png')
                plt.show()
                plt.close()

    except:
        print("Pressure coefficient not available in provided CSV file.\nRespective plot will be skipped.")

    try:
        my_cf_x = df["Skin_Friction_Coefficient_x"].tolist()
        print(
            'The maximum value of x-coordinate friction coefficient is (' + str(round(max(my_cf_x), 3)) + ') at ' + str(
                round(my_x[np.argmax(my_cf_x)], 3)) + ' m.\n')
        print(
            'The minimum value of y-coordinate friction coefficient is (' + str(round(min(my_cf_x), 3)) + ') at ' + str(
                round(my_x[np.argmin(my_cf_x)], 3)) + ' m.\n\n')
        plt.scatter(my_x, my_cf_x, c='black', s=my_size, label=my_label)
        plt.suptitle('Friction Coefficient_x (x)')
        plt.ylabel('Friction coefficient_x []')
        plt.xlabel('x [m]')
        plt.grid(True)
        if min(my_cf_x) >= 0:
            min_my_cf_x = min(my_cf_x) * 0.9
        else:
            min_my_cf_x = min(my_cf_x) * 1.1
        if max(my_cf_x) >= 0:
            max_my_cf_x = max(my_cf_x) * 1.1
        else:
            max_my_cf_x = max(my_cf_x) * 0.9

        plt.ylim((min_my_cf_x, max_my_cf_x))
        plt.xlim((min(my_x) - 0.1, max(my_x) * 1.1))
        plt.legend()
        plt.savefig('Friction Coefficient_x (x).png')
        plt.show()
        plt.close()

        if EXP_DATA_CF is not None:
            with open(EXP_DATA_CF) as f:
                exp_data = [[float(row.split()[0]), float(row.split()[1])] for row in f]
                exp_data_x = [exp_data[i][0] for i in range(0, len(exp_data))]
                exp_data_cf = [exp_data[i][1] for i in range(0, len(exp_data))]

                plt.scatter(my_x, my_cf_x, c='black', s=my_size, label=my_label)
                plt.scatter(exp_data_x, exp_data_cf, c=colour_exp_cf, s=exp_size, label=exp_label)
                plt.suptitle('Pressure Coefficient(x)')
                plt.ylabel('Pressure Coefficient []')
                plt.xlabel('x [m]')
                plt.grid(True)
                if min(min(my_cf_x), min(exp_data_cf)) >= 0:
                    min_cf = min(min(my_cf_x), min(exp_data_cf)) * 0.95
                else:
                    min_cf = min(min(my_cf_x), min(exp_data_cf)) * 1.05
                if max(max(my_cf_x), max(exp_data_cf)) >= 0:
                    max_cf = max(max(my_cf_x), max(exp_data_cf)) * 1.05
                else:
                    max_cf = max(max(my_cf_x), max(exp_data_cf)) * 0.95

                plt.ylim((my_cf_x, my_cf_x))
                plt.xlim((min(my_x) - 0.1, max(max(my_x), max(exp_data_x)) * 1.1))
                plt.legend()
                plt.savefig('Cf(x): comparison.png')
                plt.show()
                plt.close()


    except:
        print("X-coordinate Friction coefficient not available in provided CSV file.\nRespective plot will be skipped.")

    try:
        my_cf_y = df["Skin_Friction_Coefficient_y"].tolist()
        print(
            'The maximum value of y-coordinate friction coefficient is (' + str(round(max(my_cf_y), 3)) + ') at ' + str(
                round(my_x[np.argmax(my_cf_y)], 3)) + ' m.\n')
        print(
            'The minimum value of y-coordinate friction coefficient is (' + str(round(min(my_cf_y), 3)) + ') at ' + str(
                round(my_x[np.argmin(my_cf_y)], 3)) + ' m.\n\n')
        plt.scatter(my_x, my_cf_y, c='black', s=my_size, label=my_label)
        plt.suptitle('Friction Coefficient_y (x)')
        plt.ylabel('Friction coefficient_y []')
        plt.xlabel('x [m]')
        plt.grid(True)
        if min(my_cf_y) >= 0:
            min_my_cf_y = min(my_cf_y) * 0.9
        else:
            min_my_cf_y = min(my_cf_y) * 1.1
        if max(my_cf_y) >= 0:
            max_my_cf_y = max(my_cf_y) * 1.1
        else:
            max_my_cf_y = max(my_cf_y) * 0.9

        plt.ylim((min_my_cf_y, max_my_cf_y))
        plt.xlim((min(my_x) - 0.1, max(my_x) * 1.1))
        plt.savefig('Friction Coefficient_y (x).png')
        plt.show()
        plt.close()

    except:
        print("Y-coordinate Friction coefficient not available in provided CSV file.\nRespective plot will be skipped.")

    try:
        my_heat_flux = df["Heat_Flux"].tolist()
        print('The maximum value of Heat Flux is (' + str(round(max(my_heat_flux), 3)) + ') at ' + str(
            round(my_x[np.argmax(my_heat_flux)], 3)) + ' m.\n')
        print('The minimum value of Heat Flux is (' + str(round(min(my_heat_flux), 3)) + ') at ' + str(
            round(my_x[np.argmin(my_heat_flux)], 3)) + ' m.\n\n')
        plt.scatter(my_x, my_heat_flux, c='black', s=my_size, label=my_label)
        plt.suptitle('Heat Flux (x)')
        plt.ylabel('Heat Flux [$W/m^2$]')
        plt.xlabel('x [m]')
        plt.grid(True)
        if min(my_heat_flux) >= 0:
            min_my_heat_flux = min(my_heat_flux) * 0.98
        else:
            min_my_heat_flux = min(my_heat_flux) * 1.02
        if max(my_heat_flux) >= 0:
            max_my_heat_flux = max(my_heat_flux) * 1.02
        else:
            max_my_heat_flux = max(my_heat_flux) * 0.98

        plt.ylim((min_my_heat_flux, max_my_heat_flux))
        plt.xlim((min(my_x) - 0.1, max(my_x) * 1.1))
        plt.legend()
        plt.savefig('Heat Flux (x).png')
        plt.show()
        plt.close()

    except:
        print("Heat Flux not available in provided CSV file.\nRespective plot will be skipped.")

    try:
        my_y_plus = df["Y_Plus"].tolist()
        print('The maximum value of y' + get_super('+') + ' is (' + str(round(max(my_y_plus), 3)) + ') at ' + str(
            round(my_x[np.argmax(my_y_plus)], 3)) + ' m.\n')
        print('The minimum value of y' + get_super('+') + ' is (' + str(round(min(my_y_plus), 3)) + ') at ' + str(
            round(my_x[np.argmin(my_y_plus)], 3)) + ' m.\n\n')
        plt.scatter(my_x, my_y_plus, c='black', s=my_size, label=my_label)
        plt.suptitle('$y^+$(x)')
        plt.ylabel('$y^+$ []')
        plt.xlabel('x [m]')
        plt.grid(True)
        if min(my_y_plus) >= 0:
            min_my_y_plus = min(my_y_plus) * 0.9
        else:
            min_my_y_plus = min(my_y_plus) * 1.1
        if max(my_y_plus) >= 0:
            max_my_y_plus = max(my_y_plus) * 1.1
        else:
            max_my_y_plus = max(my_y_plus) * 0.9

        plt.ylim((min_my_y_plus, max_my_y_plus))
        plt.xlim((min(my_x) - 0.1, max(my_x) * 1.1))
        plt.legend()
        plt.savefig('y_plus (x).png')
        plt.show()
        plt.close()

    except:
        print("Dimensionless wall distance not available in provided CSV file.\nRespective plot will be skipped.")

    return []
