#!/usr/bin/env python3
# *************************************************************************
# Import libraries

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# *************************************************************************
def mutations_vs_angrange(df, mut_column, x_axis, directory, name):

    # .corr() returns the correlation between two columns
    pearson_a = df[mut_column].corr(df['angle_range'])

    plt.figure()

    color_values = 'burlywood'
    color_bf_line = 'rebeccapurple'

    x = df[mut_column]
    y = df['angle_range']

    m, b = np.polyfit(x, y, 1)
    plt.plot(x, m * x + b, color=color_bf_line,
             linestyle='dashed', linewidth=1)
    bf_line = 'y={:.3f}x+{:.3f}'.format(m, b)

    plt.scatter(x, y, s=3, color=color_values)

    axes = plt.gca()

    # Sets the maximum and minimum values for the axes
    # axes.autoscale(tight=True)
    axes.set_xlim([0, 60])
    axes.set_ylim([-5, 45])

    # axes.axline((0, 0), (1, 1), color='k')

    # Sets the axes labels
    plt.xlabel(f'{x_axis} mutations from germline')
    plt.ylabel('Range of the packing angle')

    # Adds graph annotations
    plt.text(s=f'Correlation: {pearson_a}', x=5, y=40, fontsize=8)
    plt.text(s=f'Best fit: {bf_line}',
             x=5, y=38, fontsize=8, color=color_bf_line)

    # plt.tight_layout()

    # Exports the figure as a .jpg file
    path_fig = os.path.join(directory, f'{name}.jpg')
    plt.savefig(path_fig, format='jpg')
    plt.close()

    return m, b
