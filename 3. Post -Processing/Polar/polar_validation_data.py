
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def polar_validation_data(AOA_old, cl_old, cd_old, cm_old, EXP_DATA, CFD_label, CFD_scatter_size, CFD_line_width, CFD_colour, exp_label, exp_line_width, exp_colour):
    if AOA_old != [None]:

        if max(AOA_old) >= 0 and min(AOA_old) >= 0:
            if min(AOA_old) >= 1:
                lim_aoa = [min(AOA_old) * 0.9, max(AOA_old)*1.1]
            else:
                lim_aoa = [-0.5, max(AOA_old) * 1.1]

        elif max(AOA_old) >= 0 and min(AOA_old) < 0:
            if max(AOA_old) >= 1:
                lim_aoa = [min(AOA_old) * 1.1, max(AOA_old) * 1.1]
            else:
                lim_aoa = [min(AOA_old) * 1.1, 1.1]
        else:
            if max(AOA_old) >= -1:
                lim_aoa = [min(AOA_old) * 1.1, 0]
            else:
                lim_aoa = [min(AOA_old) * 1.1, max(AOA_old) * 1.1]

    if cl_old != [None]:

        if max(cl_old) >= 0 and min(cl_old) >= 0:
            if min(cl_old) >= 0.1:
                lim_cl = [min(cl_old) * 0.9, max(cl_old) * 1.1]
            else:
                lim_cl = [-0.1, max(cl_old) * 1.1]

        elif max(cl_old) >= 0 and min(cl_old) < 0:
            if max(cl_old) >= 0.1 and min(cl_old) < -0.1:
                lim_cl = [min(cl_old) * 1.1, max(cl_old) * 1.1]
            elif max(cl_old) <= 0.1 and min(cl_old) < -0.1:
                lim_cl = [min(cl_old) * 1.1, 0.1]
            elif max(cl_old) >= 0.1 and min(cl_old) >= -0.1:
                lim_cl = [-0.1, max(cl_old) * 1.1]
            else:
                lim_cl = [-0.1, 0.1]

        else:
            lim_cl = [min(cl_old) * 1.1, max(cl_old) * 0.9]

    if cd_old != [None]:
        lim_cd = [0, max(cd_old) + 0.01]

    if cm_old != [None]:

        if max(cm_old) >= 0 and min(cm_old) >= 0:
            if min(cl_old) >= 0.1:
                lim_cm = [min(cm_old) * 0.9, max(cm_old) * 1.1]
            else:
                lim_cm = [-0.1, max(cm_old) * 1.1]

        elif max(cm_old) >= 0 and min(cm_old) < 0:
            if max(cm_old) >= 1:
                lim_cm = [min(cm_old) * 1.1, max(cm_old) * 1.1]
            else:
                lim_cm = [-0.1, max(cm_old) * 1.1]
        else:
            lim_cm = [min(cm_old) * 1.1, max(cm_old) * 0.9]


    if cd_old != [None] and cl_old != [None]:
        CD_CL = interpolate.interp1d(cl_old, cd_old)
        cl1 = np.arange(min(cl_old), max(cl_old), 0.001)
        cd1 = CD_CL(cl1)  # use interpolation function returned by `interp1d`

    if AOA_old != [None] and cl_old != [None]:
        CL_AOA = interpolate.interp1d(AOA_old, cl_old)
        aoa2 = np.arange(min(AOA_old), max(AOA_old), 0.01)
        cl2 = CL_AOA(aoa2)  # use interpolation function returned by `interp1d`

    if AOA_old != [None] and cm_old != [None]:
        CM_AOA = interpolate.interp1d(AOA_old, cm_old)
        aoa3 = np.arange(min(AOA_old), max(AOA_old), 0.01)
        cm3 = CM_AOA(aoa3)  # use interpolation function returned by `interp1d`

    if AOA_old != [None] and cd_old != [None]:
        CD_AOA = interpolate.interp1d(AOA_old, cd_old)
        aoa4 = np.arange(min(AOA_old), max(AOA_old), 0.01)
        cd4 = CD_AOA(aoa4)  # use interpolation function returned by `interp1d`


    if type(EXP_DATA) != str:
        print(
            'Error: no valid input for EXP_DATA. \n This item should be a string which defines the .dat file selected for reference experimental data. \n')
        sys.exit()

    elif EXP_DATA != [None]:
        with open(EXP_DATA) as f:

            if AOA_old != [None] and cd_old != [None] and cl_old != [None] and cm_old != [None]:
                exp_data = [[float(column.split(",")[0]), float(column.split(",")[1]), float(column.split(",")[2]), float(column.split(",")[3])] for column in f]
                exp_data_aoa = [exp_data[i][0] for i in range(0, len(exp_data))]
                exp_data_cd = [exp_data[i][1] for i in range(0, len(exp_data))]
                exp_data_cl = [exp_data[i][2] for i in range(0, len(exp_data))]
                exp_data_cm = [exp_data[i][3] for i in range(0, len(exp_data))]
            elif AOA_old != [None] and cd_old != [None] and cl_old != [None] and cm_old == [None]:
                exp_data = [[float(column.split(",")[0]), float(column.split(",")[1]), float(column.split(",")[2])] for column in f]
                exp_data_aoa = [exp_data[i][0] for i in range(0, len(exp_data))]
                exp_data_cd = [exp_data[i][1] for i in range(0, len(exp_data))]
                exp_data_cl = [exp_data[i][2] for i in range(0, len(exp_data))]
                exp_data_cm = None
            elif AOA_old != [None] and cd_old != [None] and cl_old == [None] and cm_old != [None]:
                exp_data = [[float(column.split(",")[0]), float(column.split(",")[1]), float(column.split(",")[2])] for column in f]
                exp_data_aoa = [exp_data[i][0] for i in range(0, len(exp_data))]
                exp_data_cd = [exp_data[i][1] for i in range(0, len(exp_data))]
                exp_data_cl = None
                exp_data_cm = [exp_data[i][2] for i in range(0, len(exp_data))]
            elif AOA_old != [None] and cd_old == [None] and cl_old != [None] and cm_old != [None]:
                exp_data = [[float(column.split(",")[0]), float(column.split(",")[1]), float(column.split(",")[2])] for column in f]
                exp_data_aoa = [exp_data[i][0] for i in range(0, len(exp_data))]
                exp_data_cd = None
                exp_data_cl = [exp_data[i][1] for i in range(0, len(exp_data))]
                exp_data_cm = [exp_data[i][2] for i in range(0, len(exp_data))]
            elif AOA_old == [None] and cd_old != [None] and cl_old != [None] and cm_old != [None]:
                exp_data = [[float(column.split(",")[0]), float(column.split(",")[1]), float(column.split(",")[2])] for column in f]
                exp_data_aoa = None
                exp_data_cd = [exp_data[i][0] for i in range(0, len(exp_data))]
                exp_data_cl = [exp_data[i][1] for i in range(0, len(exp_data))]
                exp_data_cm = [exp_data[i][2] for i in range(0, len(exp_data))]
            elif AOA_old == [None] and cd_old != [None] and cl_old != [None] and cm_old == [None]:
                exp_data = [[float(column.split(",")[0]), float(column.split(",")[1])] for column in f]
                exp_data_aoa = None
                exp_data_cd = [exp_data[i][0] for i in range(0, len(exp_data))]
                exp_data_cl = [exp_data[i][1] for i in range(0, len(exp_data))]
                exp_data_cm = None
            elif AOA_old != [None] and cd_old == [None] and cl_old != [None] and cm_old == [None]:
                exp_data = [[float(column.split(",")[0]), float(column.split(",")[1])] for column in f]
                exp_data_aoa = [exp_data[i][0] for i in range(0, len(exp_data))]
                exp_data_cd = None
                exp_data_cl = [exp_data[i][1] for i in range(0, len(exp_data))]
                exp_data_cm = None
            elif AOA_old != [None] and cd_old == [None] and cl_old == [None] and cm_old != [None]:
                exp_data = [[float(column.split(",")[0]), float(column.split(",")[1])] for column in f]
                exp_data_aoa = [exp_data[i][0] for i in range(0, len(exp_data))]
                exp_data_cd = None
                exp_data_cl = None
                exp_data_cm = [exp_data[i][1] for i in range(0, len(exp_data))]
            elif AOA_old != [None] and cd_old != [None] and cl_old == [None] and cm_old == [None]:
                exp_data = [[float(column.split(",")[0]), float(column.split(",")[1])] for column in f]
                exp_data_aoa = [exp_data[i][0] for i in range(0, len(exp_data))]
                exp_data_cd = [exp_data[i][1] for i in range(0, len(exp_data))]
                exp_data_cl = None
                exp_data_cm = None
            else:
                print('Error: option not readable by the script script')
                sys.exit()


            # Interpolation and Plot of CD(CL)
            if cd_old != [None] and cl_old != [None]:
                EXP_CD_CL = interpolate.interp1d(exp_data_cl, exp_data_cd)

                mincl = max(min(exp_data_cl), min(cl_old))
                maxcl = min(max(exp_data_cl), max(cl_old))

                CL1 = np.arange(mincl, maxcl, 0.001)

                CD1 = EXP_CD_CL(CL1)  # use interpolation function returned by `interp1d`

                err_CD_CL = [None] * len(CL1)

                for j in range(0, len(CL1)):
                    err_CD_CL[j] = 100 * (cd1[j] - CD1[j]) / CD1[j]

                if min(err_CD_CL) > 0:
                    min_err_CD_CL = 0
                else:
                    min_err_CD_CL = min(err_CD_CL) * 1.1

                if max(err_CD_CL) > 0:
                    max_err_CD_CL = max(err_CD_CL) * 1.1
                else:
                    max_err_CD_CL = 0

                plt.scatter(cl_old, cd_old, c=CFD_colour, s=CFD_scatter_size, label=CFD_label)
                plt.plot(exp_data_cl, exp_data_cd, '--', c=exp_colour, linewidth=exp_line_width, label=exp_label)
                plt.suptitle('$C_D(C_L)$: comparison with experimental data')
                plt.ylabel('Drag coefficient []')
                plt.xlabel('Lift coefficient []')
                plt.grid(True)
                plt.ylim((lim_cd[0], lim_cd[1]))
                plt.xlim((lim_cl[0], lim_cl[1]))
                plt.legend()
                plt.savefig('CD_CL_comparison1.png')
                plt.show()
                plt.close()

                plt.plot(cl1, cd1, '-', c=CFD_colour, linewidth=CFD_line_width, label=CFD_label)
                plt.plot(CL1, CD1, '--', c=exp_colour, linewidth=exp_line_width, label=exp_label)
                plt.suptitle('$C_D(C_L)$: comparison with experimental data')
                plt.ylabel('Drag coefficient []')
                plt.xlabel('Lift coefficient []')
                plt.grid(True)
                plt.ylim((lim_cd[0], lim_cd[1]))
                plt.xlim((lim_cl[0], lim_cl[1]))
                plt.legend()
                plt.savefig('CD_CL_comparison2.png')
                plt.show()
                plt.close()

                plt.scatter(CL1, err_CD_CL, c=CFD_colour, s=CFD_scatter_size, label='Error ' + CFD_label)
                plt.suptitle('$C_D(C_L)$: Error (%)')
                plt.grid(True)
                plt.ylabel('$Error_D$ [%]', size=10)
                plt.xlabel('Lift coefficient []', size=10)
                plt.ylim((min_err_CD_CL, max_err_CD_CL))
                plt.xlim((lim_cl[0], lim_cl[1]))
                plt.legend()
                plt.savefig('Error CD(CL).png')
                plt.show()
                plt.close()

            # Interpolation and Plot of CL(AOA)
            if AOA_old != [None] and cl_old != [None]:
                EXP_CL_AOA = interpolate.interp1d(exp_data_aoa, exp_data_cl)
                minaoa = max(min(exp_data_aoa), min(AOA_old))
                maxaoa = min(max(exp_data_aoa), max(AOA_old))

                AOA2 = np.arange(minaoa, maxaoa, 0.01)

                CL2 = EXP_CL_AOA(AOA2)  # use interpolation function returned by `interp1d`

                err_CL_AOA = [None] * len(AOA2)
                for i in range(0, len(AOA2)):
                    if CL2[i] == 0:
                        CL2[i] = 0.000001
                    else:
                        pass
                    err_CL_AOA[i] = 100 * (cl2[i] - CL2[i]) / CL2[i]

                if min(err_CL_AOA) > 0:
                    min_err_CL_AOA = 0
                else:
                    min_err_CL_AOA = min(err_CL_AOA) * 1.1

                if max(err_CL_AOA) > 0:
                    max_err_CL_AOA = max(err_CL_AOA) * 1.1
                else:
                    max_err_CL_AOA = 0

                plt.scatter(AOA_old, cl_old, c=CFD_colour, s=CFD_scatter_size, label=CFD_label)
                plt.plot(exp_data_aoa, exp_data_cl, '--', c=exp_colour, linewidth=exp_line_width, label=exp_label)
                plt.suptitle('$C_L$(AOA): comparison with experimental data')
                plt.ylabel('Lift coefficient []')
                plt.xlabel('Angle of attack [deg]')
                plt.grid(True)
                plt.ylim((lim_cl[0], lim_cl[1]))
                plt.xlim((lim_aoa[0], lim_aoa[1]))
                plt.legend()
                plt.savefig('CL_AOA_comparison1.png')
                plt.show()
                plt.close()

                plt.plot(aoa2, cl2, '-', c=CFD_colour, linewidth=exp_line_width, label=CFD_label)
                plt.plot(AOA2, CL2, '--', c=exp_colour, linewidth=exp_line_width, label=exp_label)
                plt.suptitle('$C_L$(AOA): comparison with experimental data')
                plt.ylabel('Lift coefficient []')
                plt.xlabel('Angle of attack [deg]')
                plt.grid(True)
                plt.ylim((lim_cl[0], lim_cl[1]))
                plt.xlim((lim_aoa[0], lim_aoa[1]))
                plt.legend()
                plt.savefig('CL_AOA_comparison2.png')
                plt.show()
                plt.close()

                plt.scatter(AOA2, err_CL_AOA, c=CFD_colour, s=CFD_scatter_size, label='Error ' + CFD_label)
                plt.suptitle('$C_L$(AOA): Error (%)')
                plt.grid(True)
                plt.ylabel('$Error_L$ [%]', size=10)
                plt.xlabel('Angle of attack [deg]', size=10)
                plt.ylim((min_err_CL_AOA, max_err_CL_AOA))
                plt.xlim((lim_aoa[0], lim_aoa[1]))
                plt.legend()
                plt.savefig('Error CL(AOA)1.png')
                plt.show()
                plt.close()

                plt.scatter(AOA2, err_CL_AOA, c=CFD_colour, s=CFD_scatter_size, label='Error ' + CFD_label)
                plt.suptitle('$C_L$(AOA): Error (%)')
                plt.grid(True)
                plt.ylabel('$Error_L$ [%]', size=10)
                plt.xlabel('Angle of attack [deg]', size=10)
                plt.ylim((-sum(err_CL_AOA) / len(err_CL_AOA)*2, sum(err_CL_AOA) / len(err_CL_AOA)*2))
                plt.xlim((lim_aoa[0], lim_aoa[1]))
                plt.legend()
                plt.savefig('Error CL(AOA)2.png')
                plt.show()
                plt.close()

            # Interpolation and Plot of CM(AOA)
            if AOA_old != [None] and cm_old != [None]:
                EXP_CM_AOA = interpolate.interp1d(exp_data_aoa, exp_data_cm)

                minaoa = max(min(exp_data_aoa), min(AOA_old))
                maxaoa = min(max(exp_data_aoa), max(AOA_old))

                AOA3 = np.arange(minaoa, maxaoa, 0.01)

                CM3 = EXP_CM_AOA(AOA3)  # use interpolation function returned by `interp1d`

                err_CM_AOA = [None] * len(AOA3)
                for i in range(0, len(AOA3)):
                    if CM3[i] == 0:
                        CM3[i] = 0.000001
                    else:
                        pass
                    err_CM_AOA[i] = 100 * (cm3[i] - CM3[i]) / CM3[i]

                if min(err_CM_AOA) > 0:
                    min_err_CM_AOA = 0
                else:
                    min_err_CM_AOA = min(err_CM_AOA) * 1.1

                if max(err_CM_AOA) > 0:
                    max_err_CM_AOA = max(err_CM_AOA) * 1.1
                else:
                    max_err_CM_AOA = 0

                plt.scatter(AOA_old, cm_old, c=CFD_colour, s=CFD_scatter_size, label=CFD_label)
                plt.plot(exp_data_aoa, exp_data_cm, '--', c=exp_colour, linewidth=exp_line_width, label=exp_label)
                plt.suptitle('$C_M$(AOA): comparison with experimental data')
                plt.ylabel('Moment coefficient []')
                plt.xlabel('Angle of attack [deg]')
                plt.grid(True)
                plt.ylim((lim_cm[0], lim_cm[1]))
                plt.xlim((lim_aoa[0], lim_aoa[1]))
                plt.legend()
                plt.savefig('CM_AOA_comparison1.png')
                plt.show()
                plt.close()

                plt.plot(aoa3, cm3, '-', c=CFD_colour, linewidth=CFD_line_width, label=CFD_label)
                plt.plot(AOA3, CM3, '--', c=exp_colour, label=exp_label)
                plt.suptitle('$C_M$(AOA): comparison with experimental data')
                plt.ylabel('Moment coefficient []')
                plt.xlabel('Angle of attack [deg]')
                plt.grid(True)
                plt.ylim((lim_cm[0], lim_cm[1]))
                plt.xlim((lim_aoa[0], lim_aoa[1]))
                plt.legend()
                plt.savefig('CM_AOA_comparison2.png')
                plt.show()
                plt.close()

                plt.scatter(AOA3, err_CM_AOA, c=CFD_colour, s=CFD_scatter_size, label='Error ' + CFD_label)
                plt.suptitle('$C_M$(AOA): Error (%)')
                plt.grid(True)
                plt.ylabel('$Error_M$ [%]', size=10)
                plt.xlabel('Angle of attack [deg]', size=10)
                plt.ylim((min_err_CM_AOA, max_err_CM_AOA))
                plt.xlim((lim_aoa[0], lim_aoa[1]))
                plt.legend()
                plt.savefig('Error CM(AOA).png')
                plt.show()
                plt.close()

            # Interpolation and Plot of CD(AOA)
            if AOA_old != [None] and cd_old != [None]:
                EXP_CD_AOA = interpolate.interp1d(exp_data_aoa, exp_data_cd)

                minaoa = max(min(exp_data_aoa), min(AOA_old))
                maxaoa = min(max(exp_data_aoa), max(AOA_old))

                AOA4 = np.arange(minaoa, maxaoa, 0.01)

                CD4 = EXP_CD_AOA(AOA4)  # use interpolation function returned by `interp1d`

                err_CD_AOA = [None] * len(AOA4)
                for i in range(0, len(AOA4)):
                    err_CD_AOA[i] = 100 * (cd4[i] - CD4[i]) / CD4[i]

                if min(err_CD_AOA) > 0:
                    min_err_CD_AOA = 0
                else:
                    min_err_CD_AOA = min(err_CD_AOA) * 1.1

                if max(err_CD_AOA) > 0:
                    max_err_CD_AOA = max(err_CD_AOA) * 1.1
                else:
                    max_err_CD_AOA = 0

                plt.scatter(AOA_old, cd_old, c=CFD_colour, s=CFD_scatter_size, label=CFD_label)
                plt.plot(exp_data_aoa, exp_data_cd, '--', c=exp_colour, linewidth=exp_line_width, label=exp_label)
                plt.suptitle('$C_D$(AOA): comparison with experimental data')
                plt.ylabel('Drag coefficient []')
                plt.xlabel('Angle of attack [deg]')
                plt.grid(True)
                plt.ylim((lim_cd[0], lim_cd[1]))
                plt.xlim((lim_aoa[0], lim_aoa[1]))
                plt.legend()
                plt.savefig('CD_AOA_comparison1.png')
                plt.show()
                plt.close()

                plt.plot(aoa4, cd4, '-', c=CFD_colour, linewidth=CFD_line_width, label=CFD_label)
                plt.plot(AOA4, CD4, '--', c=exp_colour, linewidth=exp_line_width, label=exp_label)
                plt.suptitle('$C_D$(AOA): comparison with experimental data')
                plt.ylabel('Drag coefficient []')
                plt.xlabel('Angle of attack [deg]')
                plt.grid(True)
                plt.ylim((lim_cd[0], lim_cd[1]))
                plt.xlim((lim_aoa[0], lim_aoa[1]))
                plt.legend()
                plt.savefig('CD_AOA_comparison2.png')
                plt.show()
                plt.close()

                plt.scatter(AOA4, err_CD_AOA, c=CFD_colour, s=CFD_scatter_size, label='Error ' + CFD_label)
                plt.suptitle('CD(AOA): Error (%)')
                plt.grid(True)
                plt.ylabel('$Error_D$ [%]', size=10)
                plt.xlabel('Angle of attack [deg]', size=10)
                plt.ylim((min_err_CD_AOA, max_err_CD_AOA))
                plt.xlim((lim_aoa[0], lim_aoa[1]))
                plt.legend()
                plt.savefig('Error CD(AOA).png')
                plt.show()
                plt.close()
    else:
        pass

    if cd_old != [None] and cl_old != [None]:
        plt.plot(cl1, cd1, '--', c=CFD_colour, linewidth=CFD_line_width, label=CFD_label)
        plt.suptitle('$C_D(C_L)$')
        plt.ylabel('Drag coefficient []')
        plt.xlabel('Lift coefficient []')
        plt.grid(True)
        plt.ylim((lim_cd[0], lim_cd[1]))
        plt.xlim((lim_cl[0], lim_cl[1]))
        plt.legend()
        plt.savefig('CD_CL1.png')
        plt.show()
        plt.close()

        plt.scatter(cl_old, cd_old, c=CFD_colour, s=CFD_scatter_size, label=CFD_label)
        plt.suptitle('$C_D(C_L)$')
        plt.ylabel('Drag coefficient []')
        plt.xlabel('Lift coefficient []')
        plt.grid(True)
        plt.ylim((lim_cd[0], lim_cd[1]))
        plt.xlim((lim_cl[0], lim_cl[1]))
        plt.legend()
        plt.savefig('CD_CL2.png')
        plt.show()
        plt.close()

    if AOA_old != [None] and cl_old != [None]:
        plt.plot(aoa2, cl2, '--', c=CFD_colour, linewidth=CFD_line_width, label=CFD_label)
        plt.suptitle('$C_L$(AOA)')
        plt.ylabel('Lift coefficient []')
        plt.xlabel('Angle of attack [deg]')
        plt.grid(True)
        plt.ylim((lim_cl[0], lim_cl[1]))
        plt.xlim((lim_aoa[0], lim_aoa[1]))
        plt.legend()
        plt.savefig('CL_AOA1.png')
        plt.show()
        plt.close()

        plt.scatter(AOA_old, cl_old, c=CFD_colour, s=CFD_scatter_size, label=CFD_label)
        plt.suptitle('$C_L$(AOA)')
        plt.ylabel('Lift coefficient []')
        plt.xlabel('Angle of attack [deg]')
        plt.grid(True)
        plt.ylim((lim_cl[0], lim_cl[1]))
        plt.xlim((lim_aoa[0], lim_aoa[1]))
        plt.legend()
        plt.savefig('CL_AOA2.png')
        plt.show()
        plt.close()

    if AOA_old != [None] and cm_old != [None]:
        plt.plot(aoa3, cm3, '--', c=CFD_colour, linewidth=CFD_line_width, label=CFD_label)
        plt.suptitle('$C_M$(AOA)')
        plt.ylabel('Moment coefficient []')
        plt.xlabel('Angle of attack [deg]')
        plt.grid(True)
        plt.ylim((lim_cm[0], lim_cm[1]))
        plt.xlim((lim_aoa[0], lim_aoa[1]))
        plt.legend()
        plt.savefig('CM_AOA1.png')
        plt.show()
        plt.close()

        plt.scatter(AOA_old, cm_old, c=CFD_colour, s=CFD_scatter_size, label=CFD_label)
        plt.suptitle('$C_M$(AOA)')
        plt.ylabel('Moment coefficient []')
        plt.xlabel('Angle of attack [deg]')
        plt.grid(True)
        plt.ylim((lim_cm[0], lim_cm[1]))
        plt.xlim((lim_aoa[0], lim_aoa[1]))
        plt.legend()
        plt.savefig('CM_AOA2.png')
        plt.show()
        plt.close()

    if AOA_old != [None] and cd_old != [None]:
        plt.plot(aoa4, cd4, '--', c=CFD_colour, linewidth=CFD_line_width, label=CFD_label)
        plt.suptitle('$C_D$(AOA)')
        plt.ylabel('Drag coefficient []')
        plt.xlabel('Angle of attack [deg]')
        plt.grid(True)
        plt.ylim((lim_cd[0], lim_cd[1]))
        plt.xlim((lim_aoa[0], lim_aoa[1]))
        plt.legend()
        plt.savefig('CD_AOA1.png')
        plt.show()
        plt.close()

        plt.scatter(AOA_old, cd_old, c=CFD_colour, s=CFD_scatter_size, label=CFD_label)
        plt.suptitle('$C_D$(AOA)')
        plt.ylabel('Drag coefficient []')
        plt.xlabel('Angle of attack [deg]')
        plt.grid(True)
        plt.ylim((lim_cd[0], lim_cd[1]))
        plt.xlim((lim_aoa[0], lim_aoa[1]))
        plt.legend()
        plt.savefig('CD_AOA2.png')
        plt.show()
        plt.close()

    return []
