# PARAMETRIC GENERATION OF WINGS
import matplotlib.pyplot as plt
from multiGeometry import multiGeometry
from mesh_gen import mesh_gen_EULER
from mesh_gen import mesh_gen_VISCOUS
from wall_spacing import wall_spacing
import subprocess
import time
import math as mt
from main_geo_grid import *
from input_check import input_check
from main_configuration_file import *
from lift_driver_mode import *
from cfg_printing import cfg_printing

input_check(params, alpha, dist, crel, typ, TE, Re, mu, nodes, y_plus, thick, progr, temp, ellipse_dimension, flow, Mach, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, wall_refining, ellipse_refining, ref_airfoil, n, farfield_size)
# Script General info

print('Creating airfoils geometry...')
start = time.process_time()

nairfoils = len(typ)

# Importing airfoils respective coordinates
x, y, index_Above1, index_Below1, index_Above2, index_Below2, TE_len, ref_airfoil = multiGeometry(nairfoils, n, params, alpha,
                                                                                     dist, crel, typ, TE, ref_airfoil)

plt.figure()

for i in range(0, nairfoils):
    plt.plot(x[i], y[i], 'k', linewidth=3)
    plt.axis('equal')

plt.savefig('airfoil_image.png')
plt.close()

print('Creating airfoils grid...')

if flow == 'VISCOUS':
    # Boundary Layer - structured region definition
    BL_thickness, norm_nodes, s, progr, Rex, norm_nodes_TE, U, rho = wall_spacing(Re, crel, mu, nairfoils, y_plus, progr, thick,
                                                                      temp, TE, Mach)

else:
    U = Mach * mt.sqrt(1.4 * 287 * temp)
    nu = U / Re
    rho = mu / nu

# Nodes parameters
nodes = float(nodes)

h = 0.001 / nodes  # wake

airfoil_nodes = [None] * nairfoils

for i in range(0, len(airfoil_nodes)):
    prop_list_a = [550* wall_refining[i], 250* wall_refining[i], 250* wall_refining[i], 300* wall_refining[i], 300* wall_refining[i]]  #[800, 350, 350, 420, 420]  airfoil
    airfoil_nodes[i] = [element_a * nodes for element_a in prop_list_a]

prop_list_e = [2.25* ellipse_refining, 4.5* ellipse_refining]  # ellipse
ellipse_nodes = [element_e * nodes for element_e in prop_list_e]

if flow == 'EULER':
    mesh_gen_EULER(x, y, nairfoils, index_Above1, index_Below1, index_Above2, index_Below2, alpha, dist, crel, TE, airfoil_nodes, ellipse_nodes, h, nodes, TE_len, ellipse_dimension, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, ref_airfoil, farfield_size)
else:
    mesh_gen_VISCOUS(x, y, nairfoils, index_Above1, index_Below1, index_Above2, index_Below2, alpha, dist, crel, TE, BL_thickness, norm_nodes, s, progr, airfoil_nodes, ellipse_nodes, h, nodes, norm_nodes_TE, TE_len, ellipse_dimension, Mesh_Algo, external_pts, wake_length, wake_progr, semicircle_dimension, semicircle_elem_factor, ref_airfoil, farfield_size)

data_id = open('Flow&grid_data.txt', 'w+')
data_id.write(
    'Number of airfoils: %d.\nFree-stream Mach: %.5f.\nReynolds number: %.5f.\nReference AOA: %.3f.\nReference temperature: %.2f.\nReference dynamic viscosity: %.10f.\n' % (
    nairfoils, Mach, Re, alpha[0], temp, mu))
data_id.write('type: ' + flow + '.\nnodes: %0.5f.\nellipse_dimension: %.5f.\nellipse_refining: %.6f.\n' % (nodes, ellipse_dimension, ellipse_refining))

if flow == 'VISCOUS':
    for i in range(0, nairfoils):
        data_id.write('y_plus[%d]: %.5f.\nprogr[%d]: %.5f.\nthick[%d]: %.3f.\nwall_refining[%d]: %.5f.\n' % (i, y_plus[i], i, progr[i], i, thick[i], i, wall_refining[i]))
    data_id.close()

print('Free-stream velocity is: %.3f m/s.' % U)
SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
rho_udm = ('kg/m3').translate(SUP)
print('Free-stream density is: %.6f ' % rho + rho_udm + '.\n')

print('Check the text file "Flow_data.txt" for your input flow data.\n\n')

print('Geometry and grid completed! (%.3f s)\n' % (time.process_time() - start))
subprocess.Popen("./Geo1.sh")
print('Check the file .geo before continuing.\n\n')

if CFD_run == 'YES':
    subprocess.call("./Geo2.sh")
    print("File .su2 created!\n ")
    REYNOLDS_NUMBER = Re
    MACH_NUMBER = Mach
    FREESTREAM_TEMPERATURE = temp
    MU_CONSTANT = mu
    AOA = alpha[0]
    MESH_FILENAME = 'Test.su2'
    var = [['KIND_TURB_MODEL', KIND_TURB_MODEL], ['KIND_TRANS_MODEL', KIND_TRANS_MODEL],
           ['FREESTREAM_TURBULENCEINTENSITY', FREESTREAM_TURBULENCEINTENSITY], ['RESTART_SOL', RESTART_SOL],
           ['MACH_NUMBER', MACH_NUMBER], ['AOA', AOA], ['SIDESLIP_ANGLE', SIDESLIP_ANGLE],
           ['FREESTREAM_TEMPERATURE', FREESTREAM_TEMPERATURE], ['GAS_CONSTANT', GAS_CONSTANT],
           ['FLUID_MODEL', FLUID_MODEL], ['GAMMA_VALUE', GAMMA_VALUE], ['REYNOLDS_NUMBER', REYNOLDS_NUMBER],
           ['REYNOLDS_LENGTH', REYNOLDS_LENGTH], ['VISCOSITY_MODEL', VISCOSITY_MODEL], ['MU_CONSTANT', MU_CONSTANT],
           ['MU_T_REF', MU_T_REF], ['SUTHERLAND_CONSTANT', SUTHERLAND_CONSTANT],
           ['CONDUCTIVITY_MODEL', CONDUCTIVITY_MODEL], ['PRANDTL_LAM', PRANDTL_LAM], ['PRANDTL_TURB', PRANDTL_TURB],
           ['REF_ORIGIN_MOMENT_X', REF_ORIGIN_MOMENT_X], ['REF_ORIGIN_MOMENT_Y', REF_ORIGIN_MOMENT_Y],
           ['REF_ORIGIN_MOMENT_Z', REF_ORIGIN_MOMENT_Z], ['REF_LENGTH', REF_LENGTH], ['REF_AREA', REF_AREA],
           ['CFL_NUMBER', CFL_NUMBER], ['CFL_ADAPT', CFL_ADAPT], ['CFL_ADAPT_PARAM', CFL_ADAPT_PARAM], ['ITER', ITER],
           ['CONV_NUM_METHOD_FLOW', CONV_NUM_METHOD_FLOW], ['MUSCL_FLOW', MUSCL_FLOW],
           ['SLOPE_LIMITER_FLOW', SLOPE_LIMITER_FLOW], ['VENKAT_LIMITER_COEFF', VENKAT_LIMITER_COEFF],
           ['JST_SENSOR_COEFF', JST_SENSOR_COEFF], ['CONV_FIELD', CONV_FIELD],
           ['CONV_RESIDUAL_MINVAL', CONV_RESIDUAL_MINVAL], ['CONV_STARTITER', CONV_STARTITER],
           ['CONV_CAUCHY_ELEMS', CONV_CAUCHY_ELEMS], ['CONV_CAUCHY_EPS', CONV_CAUCHY_EPS],
           ['MESH_FILENAME', MESH_FILENAME], ['SCREEN_OUTPUT', SCREEN_OUTPUT], ['HISTORY_OUTPUT', HISTORY_OUTPUT],
           ['OUTPUT_FILES', OUTPUT_FILES], ['TABULAR_FORMAT', TABULAR_FORMAT],
           ['SCREEN_WRT_FREQ_INNER', SCREEN_WRT_FREQ_INNER], ['OUTPUT_WRT_FREQ', OUTPUT_WRT_FREQ]]
    var2 = [['FIXED_CL_MODE', FIXED_CL_MODE], ['TARGET_CL', TARGET_CL], ['DCL_DALPHA', DCL_DALPHA],
            ['UPDATE_AOA_ITER_LIMIT', UPDATE_AOA_ITER_LIMIT], ['ITER_DCL_DALPHA', ITER_DCL_DALPHA],
            ['EVAL_DOF_DCX', EVAL_DOF_DCX]]
    cfg_printing(nairfoils, var, var2, flow)
    input("[ Press ENTER to run CFD simulation]\n")
    subp = subprocess.call("./CFD_run.sh")
