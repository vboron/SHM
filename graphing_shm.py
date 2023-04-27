#!/usr/bin/env python3
# *************************************************************************
# Import libraries

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# *************************************************************************
def mutations_vs_angrange(df, mut_column, x_axis, directory, name, max_val_df):
    # Plot all data
    # .corr() returns the correlation between two columns
    pearson_a = df[mut_column].corr(max_val_df['max_angle_range'])

    plt.figure()

    color_all = 'burlywood'
    color_top = 'rebeccapurple'

    x_all = df[mut_column]
    y_all = df['angle_range']

    y_max = max_val_df['max_angle_range']
    x_max = max_val_df[mut_column]

    axes = plt.gca()

    # Sets the maximum and minimum values for the axes
    # axes.autoscale(tight=True)
    axes.set_xlim([0, 50])
    axes.set_ylim([-0.5, 20])

    # axes.axline((0, 0), (1, 1), color='k')

    # Sets the axes labels
    plt.xlabel(f'{x_axis} mutations from germline')
    plt.ylabel('Range of the packing angle')

    # Adds graph annotations
    plt.text(s=f'Correlation: {pearson_a}', x=10, y=17, fontsize=8)

    # Plot highest values
    plt.scatter(x_all, y_all, s=3, color=color_all)
    plt.scatter(x_max, y_max, s=3, color=color_top)
    m, b = np.polyfit(x_max, y_max, 1)
    plt.plot(x_max, m * x_max + b, color=color_top,
            linestyle='dashed', linewidth=1)
    bf_line = 'y={:.3f}x+{:.3f}'.format(m, b)
    plt.text(s=f'Best fit: {bf_line}',
            x=10, y=15, fontsize=8, color=color_top)

    # Exports the figure as a .jpg file
    path_fig = os.path.join(directory, f'{name}.jpg')
    plt.savefig(path_fig, format='jpg')
    plt.close()

    # return m, b
    return


def hydrophobicity_histagram(x_values, labels, name):
    n_bins = 20
    color=['teal', 'peachpuff', 'dodgerblue', 'indigo', 'mediumorchid', 'lightpink']
    plt.figure()
    plt.hist(x_values, bins=n_bins, density=True, color=color, label=labels)
    plt.legend(prop={'size': 10})
    plt.xlabel(f'Mean change in hydrophobicity')
    plt.ylabel('Frequency')
    plt.savefig(f'{name}.jpg', format='jpg')

def introduced_hydrophobicity(x_values, labels, name):
    color=['black', 'silver', 'darkred', 'red', 'darkgreen', 'lawngreen', 'darkorange', 
           'bisque', 'darkmagenta', 'magenta']
    plt.figure()
    plt.hist(x_values, density=True, color='darkmagenta', label=labels)
    plt.legend(prop={'size': 10})
    plt.xlabel(f'X')
    plt.ylabel('Frequency')
    plt.savefig(f'{name}.jpg', format='jpg')


def introduced_fractional_hydrophobicity(x_values):
    groups = [*range(0, 50, 10)]
    min_max_vals = []
    for i in groups:
        i_int = int(i)
        min_max_vals.append([i_int, i_int+10])
    
    nbins = 10
    for min_mut, max_mut in min_max_vals:
        df = x_values[x_values['mut_count'].between(int(min_mut), int(max_mut))]

        def make_graph(x_col, color):
            plt.figure()
            plt.hist(df[f'fraction_{x_col}'], density=True, color=color, bins=nbins)
            plt.xlabel(f'Fraction of mutations which are {x_col}')
            plt.ylabel('Frequency')
            plt.savefig(f'{x_col}_{min_mut}_{max_mut}.jpg', format='jpg')
        
        make_graph('hydrophilic', 'turquoise')
        make_graph('hydrophobic', 'maroon')
