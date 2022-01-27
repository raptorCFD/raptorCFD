# -------------------------- CL DRIVER DEFINITION -----------------------------#
# Activate fixed lift mode (specify a CL instead of AoA, NO/YES)
FIXED_CL_MODE = 'NO'
#
# Target coefficient of lift for fixed lift mode (0.80 by default)
TARGET_CL = '0.8501'
#
# Estimation of dCL/dAlpha (0.2 per degree by default)
DCL_DALPHA = '0.2'
#
# Maximum number of iterations between AoA updates
UPDATE_AOA_ITER_LIMIT = '100'
#
# Number of iterations to evaluate dCL/dAlpha at the end of the simulation
ITER_DCL_DALPHA = '500'
#
# Evaluate dObjFunc/dCL during runtime (YES) or use the value stored in the
# direct solution file (NO).
EVAL_DOF_DCX = 'NO'
