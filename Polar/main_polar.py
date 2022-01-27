from polar_validation_data import polar_validation_data

# Insert your data
AOA = [0, 2, 4, 6, 8, 10, 12, 14, 16]
cl = [0.000090, 0.226878, 0.452728, 0.672762, 0.887041, 1.093954, 1.280726, 1.440163, 1.566406]
cd = [0.008839, 0.008998, 0.009478, 0.010262, 0.011620, 0.013439, 0.016142, 0.020469, 0.027619]
cm = [None]

CFD_label = "CFD results"
CFD_scatter_size = 8
CFD_line_width = 2
CFD_colour = 'black'

# Optional: experimental data comparison. Define file .dat or 'None'.
EXP_DATA = "Ladson (M=0.15, Re=6E6,  tripped 80 grit).dat"
exp_label = "Ladson, tripped 80 grit"
exp_line_width = 1
exp_colour = 'tab:red'

# ======================================
polar_validation_data(AOA, cl, cd, cm, EXP_DATA, CFD_label, CFD_scatter_size, CFD_line_width, CFD_colour, exp_label, exp_line_width, exp_colour)
