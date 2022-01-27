import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
import numpy as np

MY_FILE_1 = 'history.csv'
AOA = 10

df1 = pd.read_csv(MY_FILE_1)
iterations = df1["Inner_Iter"].tolist()
history_cd = df1['       "CD"       '].tolist()
history_cl = df1['       "CL"       '].tolist()
cauchy_cd = df1['   "Cauchy[CD]"   '].tolist()
cauchy_cl = df1['   "Cauchy[CL]"   '].tolist()
print(history_cd)
time = df1['    "Time(sec)"   '].tolist()
sum_time = [sum(time[:j]) for j in range(0, len(time))]
print('Total time requested was ' + str(round(sum_time[len(sum_time) - 1], 3)) + ' seconds.\n')
colour_cd = ['black', 'tab:red', 'tab:orange']
colour_cl = ['black', 'tab:blue', 'tab:cyan']

#=====================================================
exp_data_aoa = [-4.04,-2.14, 0.05, 2.05, 4.04, 6.09, 8.30, 10.12, 11.13, 12.12,13.08, 14.22, 15.26, 16.30, 17.13]
exp_data_cd = [.00871, .00800, .00809, .00816, .00823, .00885, .01050, .01201, .01239, .01332, .01503, .01625, .01900, .02218, .02560]
exp_data_cl = [-.4417, -.2385, -.0126,  .2125, .4316, .6546, .8873, 1.0707, 1.1685, 1.2605, 1.3455, 1.4365, 1.5129, 1.5739, 1.6116]

EXP_CD_AOA = interpolate.interp1d(exp_data_aoa, exp_data_cd)
EXP_CL_AOA = interpolate.interp1d(exp_data_aoa, exp_data_cl)

exp_aoa = np.arange(min(exp_data_aoa), max(exp_data_aoa), 0.001)

exp_cd_vec = EXP_CD_AOA(exp_aoa)
exp_cl_vec = EXP_CL_AOA(exp_aoa)

mindiff = 10000
for i in range(0, len(exp_cd_vec)):
    diff = abs(AOA - exp_aoa[i])
    if diff < mindiff:
        mindiff = diff
        exp_cd = exp_cd_vec[i]
        exp_cl = exp_cl_vec[i]

#=====================================================


err_cd = [None] * len(history_cd)
err_cl = [None] * len(history_cl)

err_rel_cd = [None] * len(history_cd)
err_rel_cl = [None] * len(history_cl)

for i in range(0, len(iterations)):
    err_cd[i] = 100 * (history_cd[i] - exp_cd) / exp_cd
    err_cl[i] = 100 * (history_cl[i] - exp_cl) / exp_cl
    if i != 0:
        err_rel_cd[i] = 100 * (history_cd[i - 1] - history_cd[i]) / history_cd[i - 1]
        err_rel_cl[i] = 100 * (history_cl[i - 1] - history_cl[i]) / history_cl[i - 1]
    else:
        err_rel_cd[i] = 100
        err_rel_cl[i] = 100


index_cd = 500
index_cl = 500

for j in range(len(err_rel_cd) - 1, 499, -1):
    if abs(err_rel_cd[j]) >= 0.1:
        index_cd = j
        break
for k in range(len(err_rel_cl) - 1, 499, -1):
    if abs(err_rel_cl[k]) >= 0.1:
        index_cl = k
        break

data_id = open('Convergence_data_' + str(AOA) + '.txt', 'w+')
data_id.write('Angle of attack: ' + str(AOA) + '.\n')
data_id.write('Iterations required: ' + str(iterations[len(iterations) - 1]) + '.\n')
data_id.write('Time requested: ' + str(round(sum_time[len(sum_time) - 1] , 3)) + ' seconds = ' + str(int(sum_time[len(sum_time) - 1] / 60)) + ' minutes.\n')
data_id.write('Lift coefficient goes under 0.1% of Relative Error at iteration: ' + str(index_cl) + '.\n')
data_id.write('Drag coefficient goes under 0.1% of Relative Error at iteration: ' + str(index_cd) + '.\n')
data_id.close()

fig, ax1 = plt.subplots()
plt.grid(True)
plt.suptitle('Drag Coefficient value and Cauchy vs. iterations')
ax1.set_xlabel('iteration')
ax1.plot(iterations, history_cd, color=colour_cd[1], label='$C_D$ history')
ax1.set_ylabel('Drag Coefficient', color=colour_cd[1])
ax1.set_ylim(min(history_cd[250:]), max(history_cd[250:]))
plt.legend()
ax1.tick_params(axis='y')
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(iterations, cauchy_cd, color=colour_cd[0])
ax2.set_ylabel('Cauchy Residual on $C_D$', color=colour_cd[0])  # we already handled the x-label with ax1
ax2.tick_params(axis='y')
ax2.set_ylim(min(cauchy_cd[250:]), max(cauchy_cd[250:]))
plt.xlim((250, max(iterations)))
fig.tight_layout(pad=2.5)  # otherwise the right y-label is slightly clipped
plt.savefig('History_cauchy_cd_%d.png' % AOA)
plt.show()
plt.close()



fig, ax1 = plt.subplots()
plt.suptitle('Drag Coefficient % error and Cauchy vs. iterations')
plt.grid(True)
ax1.set_xlabel('iteration')
ax1.plot(iterations, err_cd, color=colour_cd[2], label='$C_D$ Relative Error [%]')
ax1.set_ylabel('Relative Error on $C_D$', color=colour_cd[2])  # we already handled the x-label with ax1
ax1.tick_params(axis='y')
ax1.set_ylim(min(err_cd[250:]), max(err_cd[250:]))
plt.legend()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(iterations, cauchy_cd, color=colour_cd[0])
ax2.set_ylabel('Cauchy Residual on $C_D$', color=colour_cd[0])  # we already handled the x-label with ax1
ax2.tick_params(axis='y')
ax2.set_ylim(min(cauchy_cd[250:]), max(cauchy_cd[250:]))
plt.xlim((250, max(iterations)))
fig.tight_layout(pad=2.5)  # otherwise the right y-label is slightly cdipped
plt.savefig('Error_cauchy_cd_%d.png' % AOA)
plt.show()
plt.close()


fig, ax1 = plt.subplots()
plt.suptitle('Drag Coefficient % Relative Error vs. iterations')
plt.grid(True)
ax1.set_xlabel('iteration')
ax1.plot(iterations, err_rel_cd, color=colour_cd[2], label='$C_D$ Relative Error [%]')
ax1.set_ylabel('Relative Error on $C_D$', color=colour_cd[2])  # we already handled the x-label with ax1
ax1.tick_params(axis='y')
ax1.set_ylim(min(err_rel_cd[250:]), max(err_rel_cd[250:]))
plt.xlim((500, max(iterations)))
fig.tight_layout(pad=2.5)  # otherwise the right y-label is slightly dipped
plt.savefig('RelError_cauchy_cd_%d.png' % AOA)
plt.legend()
plt.show()
plt.close()


fig, ax1 = plt.subplots()
plt.grid(True)
plt.suptitle('Lift Coefficient value and Cauchy vs. iterations')
ax1.set_xlabel('iteration')
ax1.plot(iterations, history_cl, color=colour_cl[1], label='$C_L$ history')
ax1.set_ylabel('Lift Coefficient', color=colour_cl[1])
ax1.set_ylim(min(history_cl[250:]), max(history_cl[250:]))
ax1.tick_params(axis='y')
plt.legend()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(iterations, cauchy_cl, color=colour_cl[0])
ax2.set_ylabel('Cauchy Residual on $C_L$', color=colour_cl[0])  # we already handled the x-label with ax1
ax2.tick_params(axis='y')
ax2.set_ylim(min(cauchy_cl[250:]), max(cauchy_cl[250:]))
plt.xlim((500, max(iterations)))
fig.tight_layout(pad=2.5)  # otherwise the right y-label is slightly clipped
plt.savefig('History_cauchy_cl_%d.png' % AOA)
plt.show()
plt.close()

fig, ax1 = plt.subplots()
plt.suptitle('Lift Coefficient % error and Cauchy vs. iterations')
plt.grid(True)
ax1.set_xlabel('iteration')
ax1.plot(iterations, err_cl, color=colour_cl[2], label='$C_L$ Relative Error [%]')
ax1.set_ylabel('Error on $C_L$', color=colour_cl[2])  # we already handled the x-label with ax1
ax1.tick_params(axis='y')
ax1.set_ylim(min(err_cl[250:]), max(err_cl[250:]))
plt.legend()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(iterations, cauchy_cl, color=colour_cl[0])
ax2.set_ylabel('Cauchy Residual on $C_L$', color=colour_cl[0])  # we already handled the x-label with ax1
ax2.tick_params(axis='y')
ax2.set_ylim(min(cauchy_cl[250:]), max(cauchy_cl[250:]))
plt.xlim((500, max(iterations)))
fig.tight_layout(pad=2.5)  # otherwise the right y-label is slightly clipped
plt.savefig('Error_cauchy_cl_%d.png' % AOA)
plt.show()
plt.close()

fig, ax1 = plt.subplots()
plt.suptitle('Lift Coefficient % Relative Error vs. iterations')
plt.grid(True)
ax1.set_xlabel('iteration')
ax1.plot(iterations, err_rel_cl, color=colour_cl[2], label='$C_L$ Relative Error [%]')
ax1.set_ylabel('Relative Error on $C_L$', color=colour_cl[2])  # we already handled the x-label with ax1
ax1.tick_params(axis='y')
ax1.set_ylim(min(err_rel_cl[250:]), max(err_rel_cl[250:]))
plt.xlim((500, max(iterations)))
fig.tight_layout(pad=2.5)  # otherwise the right y-label is slightly clipped
plt.legend()
plt.savefig('RelError_cauchy_cl_%d.png' % AOA)
plt.show()
plt.close()
