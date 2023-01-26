#!/usr/bin/env python3
# *************************************************************************
# Import libraries 

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



# *************************************************************************
def mutations_vs_angrange(df, directory, name):


    # .corr() returns the correlation between two columns
    pearson_a = df['total_mut'].corr(df['angle_range'])

    plt.figure()

    color_values = 'mediumpurple'
    color_bf_line = 'rebeccapurple'

    x = df['total_mut']
    y = df['angle_range']

    m, b = np.polyfit(x, y, 1)
    plt.plot(x, m * x + b, color=color_bf_line, linestyle='dashed', linewidth=1)

    plt.scatter(x, y, s=2, color=color_values)

    axes = plt.gca()

    # Sets the maximum and minimum values for the axes
    # axes.autoscale(tight=True)
    axes.set_xlim([0, 60])
    axes.set_ylim([0, 45])

    # axes.axline((0, 0), (1, 1), color='k')

    # Sets the axes labels
    plt.xlabel('VH + VL mutations from germline')
    plt.ylabel('Range of the packing angle')

    # Adds graph annotations
    # plt.text(s='Line: y=x', x=-66, y=-27, fontsize=8)

    # plt.text(s='Best fit: y={:.3f}x+{:.3f}'.format(m, b),
    #         x=-66, y=-29, fontsize=8, color=color_bf_line)

    # plt.tight_layout()

    # add best fit data to dataframe and export the dataframe
    # add best fit lines to statistics dataframe
    bf_line = 'y={:.3f}x+{:.3f}'.format(m, b)

    # Exports the figure as a .jpg file
    path_fig = os.path.join(directory, f'{name}.jpg')
    plt.savefig(path_fig, format='jpg')
    plt.close()

    return m, b
