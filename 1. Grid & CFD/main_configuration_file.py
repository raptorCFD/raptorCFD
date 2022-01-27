#==================================================================================================================
# ------------- DIRECT, ADJOINT, AND LINEARIZED PROBLEM DEFINITION ------------#
# Specify turbulent model (NONE, SA, SA_NEG, SST)
KIND_TURB_MODEL = 'SA'

# Specify transition model (NONE, BC)
KIND_TRANS_MODEL = 'NONE'

FREESTREAM_TURBULENCEINTENSITY = '0.0015'

# Restart solution (NO, YES)
RESTART_SOL = 'NO'

# -------------------- COMPRESSIBLE FREE-STREAM DEFINITION --------------------#
# Side-slip angle (degrees, only for compressible flows)
SIDESLIP_ANGLE = '0.0'

# Reynolds length (1 m by default)
REYNOLDS_LENGTH = '1.04167'
#REYNOLDS_LENGTH = '1.3'

# ---- IDEAL GAS, POLYTROPIC, VAN DER WAALS AND PENG ROBINSON CONSTANTS -------#
# Different gas model (STANDARD_AIR, IDEAL_GAS, VW_GAS, PR_GAS)
FLUID_MODEL = 'STANDARD_AIR'

# Specific gas constant (287.058 J/kg*K default and this value is hardcoded
#                        for the model STANDARD_AIR, compressible only)
GAS_CONSTANT = '287.058'

# Ratio of specific heats (1.4 default and the value is hardcoded
#                          for the model STANDARD_AIR, compressible only)
GAMMA_VALUE = '1.4'

# --------------------------- VISCOSITY MODEL ---------------------------------#
# Viscosity model (SUTHERLAND, CONSTANT_VISCOSITY).
VISCOSITY_MODEL = 'CONSTANT_VISCOSITY'

# Sutherland Temperature Ref (273.15 K default value for AIR SI)
# (only for VISCOSITY_MODEL = SUTHERLAND)
MU_T_REF = '273.15'

# Sutherland constant (110.4 default value for AIR SI)
# (only for VISCOSITY_MODEL = SUTHERLAND)
SUTHERLAND_CONSTANT = '110.4'

# --------------------------- THERMAL CONDUCTIVITY MODEL ----------------------#
# Conductivity model (CONSTANT_CONDUCTIVITY, CONSTANT_PRANDTL).
CONDUCTIVITY_MODEL = 'CONSTANT_PRANDTL'
#
# Laminar Prandtl number (0.72 (air), only for CONSTANT_PRANDTL)
PRANDTL_LAM = 0.72
#
# Turbulent Prandtl number (0.9 (air), only for CONSTANT_PRANDTL)
PRANDTL_TURB = 0.90

# ---------------------- REFERENCE VALUE DEFINITION ---------------------------#
# Reference origin for moment computation
REF_ORIGIN_MOMENT_X = '0.25'
REF_ORIGIN_MOMENT_Y = '0.00'
REF_ORIGIN_MOMENT_Z = '0.00'
#
# Reference length for pitching, rolling, and yawing non-dimensional moment
REF_LENGTH = '1.0'

# Reference area for force coefficients (0 implies automatic calculation)
REF_AREA = '0'

# ------------- COMMON PARAMETERS DEFINING THE NUMERICAL METHOD ---------------#
# Courant-Friedrichs-Lewy condition of the finest grid
CFL_NUMBER = '500.0'

# Adaptive CFL number (NO, YES)
CFL_ADAPT = 'NO'

# Parameters of the adaptive CFL number (factor down, factor up, CFL min value,
#                                        CFL max value )
CFL_ADAPT_PARAM = '( 0.999, 100 , 100.0, 1e5 )'

# Number of total iterations
ITER = '600'

# -------------------- FLOW NUMERICAL METHOD DEFINITION -----------------------#
# Convective numerical method (JST, LAX-FRIEDRICH, CUSP, ROE, AUSM, HLLC,
#                              TURKEL_PREC, MSW)
CONV_NUM_METHOD_FLOW = 'ROE'

# Monotonic Upwind Scheme for Conservation Laws (TVD) in the flow equations.
#           Required for 2nd order upwind schemes (NO, YES)
MUSCL_FLOW = 'YES'

# Slope limiter (VENKATAKRISHNAN, MINMOD)
SLOPE_LIMITER_FLOW = 'VENKATAKRISHNAN'

# Coefficient for the Venkat's limiter (upwind scheme). A larger values decrease
#             the extent of limiting, values approaching zero cause
#             lower-order approximation to the solution (0.05 by default)
VENKAT_LIMITER_COEFF = '0.05'

# 2nd and 4th order artificial dissipation coefficients for JST method ( 0.5, 0.02 by default )
JST_SENSOR_COEFF = '( 0.5, 0.02 )'

# --------------------------- CONVERGENCE PARAMETERS --------------------------#
# Convergence criteria (CAUCHY, RESIDUAL)
CONV_FIELD= '(LIFT, DRAG, RMS_DENSITY)'
#
# Min value of the residual (log10 of the residual)
CONV_RESIDUAL_MINVAL= '-4'
#
# Start convergence criteria at iteration number
CONV_STARTITER= '500'
#
# Number of elements to apply the criteria
CONV_CAUCHY_ELEMS= '200'
#
# Epsilon to control the series convergence
CONV_CAUCHY_EPS= '5e-3'

# ------------------------- INPUT/OUTPUT INFORMATION --------------------------#
# Screen output fields
SCREEN_OUTPUT = '(INNER_ITER, WALL_TIME, RMS_DENSITY, RMS_NU_TILDE, LIFT, DRAG, AVG_CFL)'

# History output groups
HISTORY_OUTPUT= '(ITER, WALL_TIME, RMS_RES, AERO_COEFF, CAUCHY)'

# Files to output
# Possible formats : (TECPLOT, TECPLOT_BINARY, SURFACE_TECPLOT,
#  SURFACE_TECPLOT_BINARY, CSV, SURFACE_CSV, PARAVIEW, PARAVIEW_BINARY, SURFACE_PARAVIEW,
#  SURFACE_PARAVIEW_BINARY, MESH, RESTART_BINARY, RESTART_ASCII, CGNS, STL)
# default : (RESTART, PARAVIEW, SURFACE_PARAVIEW)
OUTPUT_FILES = '(RESTART, PARAVIEW, SURFACE_PARAVIEW, SURFACE_CSV)'

# Output file format (PARAVIEW, TECPLOT, STL)
TABULAR_FORMAT = 'CSV'

# Writing frequency for screen output
SCREEN_WRT_FREQ_INNER = '50'

# Writing frequency for volume/surface output
OUTPUT_WRT_FREQ = '100'


