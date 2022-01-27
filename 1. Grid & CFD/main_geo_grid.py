# ============================== TYPOLOGY
# 'YES': grid generation + CFD simulation
# 'NO': grid generation only
CFD_run = 'NO'

# ============================== INPUT DATA
# The length of each item must be the same and equal to the number of airfoils.
# Define chosen parametrization of each airfoil.
typ = ["IGP", "IGP"] #, "MY_FILE_3"]

# Specify geometric parameters.
params = [[0.68789098, 0.46936949, 0.02463844, 0.03854695, 0.30437879, 0.17412334, 1.23036161, 1.58036565], [0.27813020, 0.43055506, 0.08076999, 0.10462659, 0.35492561, 0.11903372, 0.95241605, 3.11448588]]
 
# Specify Trailing Edge type
#TE = ['N', 'Y', 'N'] #C-11
TE = ['N', 'Y']

# Absolute angle of attack of each airfoil.
#alpha = [9.081, 30.481, 39.511]
alpha = [7.7, 37.7]

# Relative gap between TE and LE of each airfoil w.r.t. previous [x, y coordinates].
#dist = [[0, 0], [-0.04396, -0.054945], [-0.03846, -0.04375]]
dist = [[0, 0], [0.0, -0.02]]

# Relative respective chord (%).
#crel = [1, 0.2442, 0.2093] #C-11
crel = [1, 0.3]

# ============================== GRID DATA
farfield_size = 300
# Grid: nodes' parameter.
nodes = 0.95

# Ellipse dimension.
ellipse_dimension = 0.85
ellipse_refining = 0.9

# Structured region
y_plus = [0.08, 0.12] #, 0.02]
thick = [0.15, 0.25] #, 1, 1]
progr = [1.1, 1.1] #, 1.1]
wall_refining = [1.25, 1.25]#, 1.25]

flow = 'VISCOUS'

#cc-1
Mach = 0.075
Re = 1e6
mu = 1.853e-05
temp = 300


# ============================== ADVANCED PARAMETERS
# Reference airfoil
# In case of multi-element configuration, define which is the airfoil to take as main reference for AOA, relative chord and grid features.
# If single-element configuration, insert 1. "DEFAULT" takes the airfoil with higher chord as reference.
ref_airfoil = "DEFAULT"

# Wake length and elements progression inside.
wake_length = 50
wake_progr = 1.3

# Meshing algorithm [an integer from 1 to 9 or 'DEFAULT', which corresponds to 5 (Delaunay)].
Mesh_Algo = 6

# Optional additional points, to help meshing algorithm (strictly depending on geometry and size).
# Note: almost necessary for Frontal-Delaunay (Mesh_Algo = 6); could avoid issues for Delaunay (Mesh_Algo = 5).
external_pts = 'YES'
# Circle dimension and elements' size (only if external_pts = 'YES')
# Dimension: value multiplied by greater ellipse's semi-axis. This product has a maximum of 200. Suggested: 15 - 25.
semicircle_dimension = 18
# Elements' size: value multiplied by minimum elements' size of the grid. Suggested: 200 - 300.
semicircle_elem_factor = 250

# Number of each airfoil's generating points. Suggested: 500.
n = 300
