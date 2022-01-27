import matplotlib.pyplot as plt
from validation_data import validation_data
import pandas as pd
import numpy as np

# CFD Results
MY_FILE = "surface_flow.csv"
my_label = "CFD Results"
my_size = 2


# Experimental Results
exp_label = "Experimental Results"
exp_size = 10

# Experimental Data (Cp)
EXP_DATA_CP = None
colour_exp_cp = 'tab:orange'

# Experimental Data (Cf)
EXP_DATA_CF = None
colour_exp_cf = 'tab:green'


# =======================================
validation_data(MY_FILE, EXP_DATA_CP, EXP_DATA_CF, my_label, exp_label, my_size, exp_size, colour_exp_cp, colour_exp_cf)
