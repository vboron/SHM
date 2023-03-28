#!/usr/bin/env python3

import sys

def check_equal(computed, expected):
    test_name = sys._getframe().f_back.f_code.co_name
    if computed == expected:
        print(f'Check within {test_name} PASSED')
    else:
        print(f'Check within {test_name} FAILED; computed={computed}, expected={expected}')
    assert computed == expected

def one_letter_code(pdb, res):
    """
    Go from the three-letter code to the one-letter code.

    Input:  residue      --- Three-letter residue identifier
    Return: one_letter   --- The one-letter residue identifier

    20.10.2020  Original   By: LD
    """

    dic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F',
           'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E',
           'TYR': 'Y', 'MET': 'M', 'XAA': 'X', 'UNK': 'X'}
    if res not in dic:
        raise ValueError("{}: {} not in dic".format(pdb, res))

    one_letter = dic[res]
    return one_letter


def hydrophobicity(resi):
    # 3. eisenberg consensus hydrophobicity
    # Consensus values: Eisenberg, et al 'Faraday Symp.Chem.Soc'17(1982)109
    Hydrophathy_index = {'A': 0.25, 'R': -1.76, "N": -0.64, "D": -0.72, "C": 0.04, "Q": -0.69, "E": -0.62,
                         "G": 0.16, "H": -0.40, "I": 0.73, "L": 0.53, "K": -1.10, "M": 0.26, "F": 0.61,
                         "P": -0.07,
                         "S": -0.26, "T": -0.18, "W": 0.37, "Y": 0.02, "V": 0.54, "X": -0.5}  # -0.5 is average
    if resi in list(Hydrophathy_index.keys()):
        hydrophobicity = Hydrophathy_index[resi]
    else:
        hydrophobicity = 0
    return hydrophobicity