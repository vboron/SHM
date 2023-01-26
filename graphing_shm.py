# *************************************************************************
# Import libraries 

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import utils


# *************************************************************************
def actual_vs_predicted_from_df(df, directory):


    # .corr() returns the correlation between two columns
    pearson_a = df['angle'].corr(df['predicted'])

    plt.figure()

    color_values = 'mediumpurple'
    color_bf_line = 'rebeccapurple'

    x = df['angle']
    y = df['predicted']

    m, b = np.polyfit(x, y, 1)
    plt.plot(x, m * x + b, color=color_bf_line, linestyle='dashed', linewidth=1)

    plt.scatter(x, y, s=2, color=color_values)

    axes = plt.gca()

    # Sets the maximum and minimum values for the axes
    # axes.autoscale(tight=True)
    axes.set_xlim([-70, -20])
    axes.set_ylim([-70, -20])

    axes.axline((0, 0), (1, 1), color='k')

    # Sets the axes labels
    plt.xlabel('Actual interface angle')
    plt.ylabel('Predicted interface angle')

    # Adds graph annotations
    plt.text(s='Line: y=x', x=-66, y=-27, fontsize=8)

    plt.text(s='Best fit: y={:.3f}x+{:.3f}'.format(m, b),
            x=-66, y=-29, fontsize=8, color=color_bf_line)
    plt.text(s='RELRMSE: {:.3}'.format(
            float(stats['RELRMSE'])), x=-66, y=-31, fontsize=8)

    plt.tight_layout()

    # add best fit data to dataframe and export the dataframe
    # add best fit lines to statistics dataframe
    bf_line = 'y={:.3f}x+{:.3f}'.format(m, b)

    stats['fit'] = bf_line
    stats['slope'] = m
    stats['intercept'] = b

    path_stats = os.path.join(directory, f'{stats_csv_name}_stats.csv')
    stats.to_csv(path_stats, index=False)

    # Exports the figure as a .jpg file
    path_fig = os.path.join(directory, f'{pa_graph_name}.jpg')
    plt.savefig(path_fig, format='jpg')
    plt.close()

    return m, b
