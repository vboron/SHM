#!/usr/bin/env python3
# *************************************************************************
# Import libraries

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
import seaborn as sns


# *************************************************************************
def mutations_vs_angrange(df, mut_column, x_axis, directory, name, max_val_df):
    # Plot all data
    # .corr() returns the correlation between two columns
    pearson_max = df[mut_column].corr(max_val_df['max_angle_range'])
    pearson_all = df[mut_column].corr(df['angle_range'])

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
    # plt.text(s=f'Correlation for max: {pearson_a:.3f}', x=10, y=17, fontsize=8)

    # Plot highest values
    plt.scatter(x_all, y_all, s=3, color=color_all)
    plt.scatter(x_max, y_max, s=3, color=color_top)

    m, b = np.polyfit(x_all, y_all, 1)
    plt.plot(x_all, m * x_all + b, color=color_all,
            linestyle='dashed', linewidth=1)
    bf_line = 'y={:.3f}x+{:.3f}'.format(m, b)
    plt.text(s=f'Best fit: {bf_line}',
            x=10, y=17, fontsize=8, color=color_all)
    m_max, b_max = np.polyfit(x_max, y_max, 1)
    plt.plot(x_max, m_max * x_max + b_max, color=color_top,
            linestyle='dashed', linewidth=1)
    bf_line_max = 'y={:.3f}x+{:.3f}'.format(m_max, b_max)
    plt.text(s=f'Best fit for max values: {bf_line_max}',
            x=10, y=15, fontsize=8, color=color_top)

    # Exports the figure as a .jpg file
    path_fig = os.path.join(directory, f'agl_{name}_{x_axis}_graph.jpg')
    plt.savefig(path_fig, format='jpg')
    plt.close()

    print(f'Pearson_all_{x_axis}: {pearson_all} \n Pearson_max_{x_axis}: {pearson_max}')
    # return m, b
    return pearson_all, pearson_max


def hydrophobicity_histagram(x_values, labels, name):
    n_bins = 20
    color=['teal', 'peachpuff', 'dodgerblue', 'indigo', 'mediumorchid', 'lightpink']
    plt.figure()
    plt.hist(x_values, bins=n_bins, density=True, color=color, label=labels)
    plt.legend(prop={'size': 10})
    plt.xlabel(f'Mean change in hydrophobicity')
    plt.ylabel('Frequency')
    plt.savefig(f'{name}.jpg', format='jpg')

def introduced_hydrophobicity(df):
    print(df)
    x = df['length_CDRs']
    y_values = df.drop(['length_CDRs'], axis=1)

    for col in y_values:
        plt.figure()
        color = ''
        label_y = ''
        if 'hydrophilic' in col:
            color = 'turquoise'
            label_y = 'Number of hydrophilic residues'
        if 'hydrophobic' in col:
            color = 'maroon'
            label_y = 'Number of hydrophobic residues'
        y = df[col]
        pearson_a = x.corr(y)
        plt.scatter(x=x, y=y, s=3, color=color)
        plt.xlabel(f'Number of mutations from germline')
        plt.ylabel(label_y)
        m, b = np.polyfit(x, y, 1)
        plt.plot(x, m * x + b, color='black',
                linestyle='dashed', linewidth=1)
        bf_line = 'y={:.3f}x+{:.3f}'.format(m, b)
        plt.text(s=f'Best fit: {bf_line}',
            x=5, y=8, fontsize=8, color='black')
        plt.text(s=f'Correlation: {pearson_a}', x=5, y=6, fontsize=8)
        axes = plt.gca()
        axes.set_xlim([0, 50])
        # Fix axes as integer
        axes.xaxis.set_major_locator(MaxNLocator(integer=True))
        axes.set_ylim([0, 10])
        axes.yaxis.set_major_locator(MaxNLocator(integer=True))
        plt.savefig(f'{col}.jpg', format='jpg')
        plt.close()


def introduced_fractional_hydrophobicity(x_values):
    groups = [*range(1, 50, 10)]
    min_max_vals = []
    for i in groups:
        i_int = int(i)
        min_max_vals.append([i_int, i_int+10])
    
    for min_mut, max_mut in min_max_vals:
        df = x_values[x_values['mut_count'].between(int(min_mut), int(max_mut))]

        def make_graph(x_col, color):
            plt.figure()
            plt.hist(df[f'fraction_{x_col}'], color=color)
            plt.xlabel(f'Fraction of mutations which are {x_col}')
            plt.ylabel('Frequency')
            axes = plt.gca()
            axes.set_xlim([0, 1])
            axes.yaxis.set_major_locator(MaxNLocator(integer=True))
            plt.savefig(f'{x_col}_{min_mut}_{max_mut}.jpg', format='jpg')
            plt.close()
        
        make_graph('hydrophilic', 'turquoise')
        make_graph('hydrophobic', 'maroon')


def hydrophobicity_change_vs_mutation(df):
    """
        cdr_dY  cdr_dH  cdr_len  cdr_mut  ...  fv_dY  fv_dH  fv_len  fv_mut
    0       1    0.88       39        9  ...      0   2.22     194      25
    1      -5   -0.10       38       21  ...     -6   1.31     196      54
    2       0    1.08       39       17  ...      1   0.79     193      26
    3       0   -0.28       39        2  ...      0  -0.12     193       4
    4      -1   -2.97       39       14  ...     -1  -6.04     194      29
    """

    metrics = []

    def make_graph(x, y, title, y_axis):
        pearson_r = x.corr(y)
        metrics.append(pearson_r)
        ax = sns.regplot(y=y, x=x, color='darkturquoise')
        ax.set(xlabel='mutations from germline', ylabel=y_axis)
        plt.savefig(f'{title}.png')
        plt.clf()
        print(f'Correlation={pearson_r}')

    for region in ['cdr', 'fwk', 'fv']:
        print(f'Graphing region {region} mutations against dY...')
        make_graph(df[f'{region}_mut'], df[f'{region}_dY'], f'{region}_dY_graph', 'change in the number of Tyrosines')
        print(f'Graphing region {region} mutations against dH...')
        make_graph(df[f'{region}_mut'], df[f'{region}_dH'], f'{region}_dH_graph', 'change in the hydrophobicity')
    
    pd.DataFrame(columns=['cdr_dY', 'cdr_dH', 'fwk_dY', 'fwk_dH', 'fv_dY', 'fv_dH'], data=[metrics]).to_csv('dYdH_pearsons.csv', index=False)
