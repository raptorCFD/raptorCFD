import numpy as np
import math as mt
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt


def geo_converter(input_filename, output_filename, rotation, scaling):
    with open(input_filename) as f:
        data = [[float(row.split()[0]), float(row.split()[1])] for row in f]

    x = np.asarray([data[i][0] for i in range(0, len(data))])
    y = np.asarray([data[i][1] for i in range(0, len(data))])

    # Input Plot
    plt.plot(x, y, 'k', linewidth=3)
    plt.axis('equal')
    plt.grid(True)
    plt.savefig('%s_input.png' % input_filename)
    plt.show()
    plt.close()


    # AOA calculus: linear regression
    x_reshape = np.array(x).reshape((-1, 1))
    y_reshape = np.array(y).reshape((-1, 1))
    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(x_reshape, y_reshape)  # perform linear regression
    y_pred = linear_regressor.predict(x_reshape)  # make predictions
    plt.scatter(x_reshape, y_reshape)
    plt.plot(x_reshape, y_pred, color='red')
    plt.show()
    plt.close()
    alpha_input = mt.atan((y_pred[int(len(y_pred)/2) - 1] - y_pred[0]) / (x[int(len(x)/2) - 1] - x[0]))

    if rotation == 'YES':
        # 1st rotation
        for i in range(0, len(x)):
            x[i] = mt.cos(-alpha_input) * x[i] - mt.sin(-alpha_input) * y[i]
            y[i] = mt.sin(-alpha_input) * x[i] + mt.cos(-alpha_input) * y[i]
    else:
        pass

    # Expansion / compression
    if scaling == 'YES':
        c_input = abs(max(x) - min(x))
        for j in range(0, len(x)):
            x[j] = x[j] / c_input
            y[j] = y[j] / c_input

    x_min = min(x)
    y_min = y[np.argmin(x)]

    # Centering on origin
    for k in range(0, len(x)):
        x[k] = x[k] - x_min
        y[k] = y[k] - y_min

    if rotation == 'YES':
        # 2nd rotation: TE coincides with x = 1, y = 0
        new_alpha = mt.atan((y[np.argmax(x)] - y[np.argmin(x)])/(x[np.argmax(x)] - x[np.argmin(x)]))
        for l in range(0, len(x)):
            x[l] = mt.cos(-new_alpha) * x[l] - mt.sin(-new_alpha) * y[l]
            y[l] = mt.sin(-new_alpha) * x[l] + mt.cos(-new_alpha) * y[l]
    else:
        pass
    plt.plot(x, y, 'k', linewidth=3)
    plt.axis('equal')
    plt.grid(True)
    plt.savefig('%s_output.png' % output_filename)
    plt.show()
    plt.close()

    # Writing text file
    data = np.array([x, y])
    data = data.T
    # here you transpose your data, so to have it in two columns

    with open(output_filename, 'w+') as datafile_id:
        # here you open the ascii file

        np.savetxt(datafile_id, data, fmt=['%.10f', '%.10f'])
        # here the ascii file is written.


    print('Conversion completed!\nCheck the file "%s" for your airfoil coordinates converted.' % (output_filename))

    return []
