#!/usr/bin/env python3

# OS   Homo sapiens (human)
# /translation=
# OS   Mus musculus (house mouse)
# AC
# FT                   /organism="Mus musculus"
# FT                   /organism="Homo sapiens"
# /protein_id=
import argparse
import os
import pandas as pd
import subprocess
import re
import graphing_shm as graph


def get_data4fasta(file, dire):
    rel_lines = ''
    with open(os.path.join(file), 'r') as f:
        for line in f:
            if line.startswith('FT                   ') or line.startswith('AC') or line.startswith('//'):
                rel_lines = rel_lines + line
    files = rel_lines.split('//')
    files = [file for file in files if 'FT                   /organism="Homo sapiens"' in file or 'FT                   /organism="Mus musculus"' in file]
    for file in files:
        filter_file(file, dire)
    return


def filter_file(f, direct):
    fasta_data = ['AC', 'FT                   /protein_id', 'FT                   /translation']
    ignore = ['FT                   /', 'AC']
    lines = f.split('\n')
    firsthalf = [l for l in lines if l.startswith(tuple(fasta_data))]
    secondhalf = [l for l in lines if not l.startswith(tuple(ignore))]
    rel_lines = firsthalf + secondhalf
    rel_lines = [l.replace(' ', '') for l in rel_lines]
    rel_lines = [l[2:].replace('"', '') for l in rel_lines if l != '']
    rel_lines = [l.replace('/protein_id=', '') for l in rel_lines]
    rel_lines = [l.replace('/translation=', '') for l in rel_lines]
    rel_lines = [l.replace(';', '') for l in rel_lines]
    rel_lines = [l for l in rel_lines if not l.islower()]
    rel_lines = [l for l in rel_lines if l.isupper() == True]
    final = rel_lines[:2] + [''.join(rel_lines[2:])]

    with open(os.path.join(direct, f'{final[0]}.faa'), 'w') as fasta:
        fasta.write(f'>{final[1]}\n{final[2]}')
    return


def run_test_get_data4fasta():
    test_input = """
    ID   JA013192; SV 1; linear; unassigned DNA; PAT; HUM; 840 BP.
XX
AC   JA013192;
XX
DT   05-APR-2011 (Rel. 108, Created)
DT   05-APR-2011 (Rel. 108, Last updated, Version 1)
XX
DE   Sequence 20 from Patent EP2278020.
XX
KW   .
XX
OS   Homo sapiens (human)
OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia;
OC   Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae;
OC   Homo.
XX
RN   [1]
RA   CLARK K., JOHNSON P.J.;
RT   "Antibody gene transfer and recombinant AAV therefor";
RL   Patent number EP2278020-A2/20, 26-JAN-2011.
RL   NATIONWIDE CHILDREN S HOSPITAL INC [US].
XX
DR   MD5; aef8177467647257a7607dea492d55b7.
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..840
FT                   /organism="Homo sapiens"
FT                   /mol_type="unassigned DNA"
FT                   /db_xref="taxon:9606"
FT   CDS             1..840
FT                   /transl_table=1
FT                   /protein_id="CCA61653.1"
FT                   /translation="MWWRLWWLLLLLLLLWPMVWADIVLTQSPGTLSLSAGERATLSCR
FT                   ASQSVSSGSLAWYQQKPGQAPRLLIYGASTRATGIPDRFSGSGSGTDFTLTIGRLEPED
FT                   LAVYYCQQYGTSPYTFGQGTKVDIKRGGGGSGGGGSGGGGSRSSQVQLVQSGAEVKKPG
FT                   SSVQVSCKASGGTFSMYGFNWVRQAPGHGLEWMGGIIPIFGTSNYAQKFRGRVTFTADQ
FT                   ATSTAYMELTNLRSDDTAVYYCARDFGPDWEDGDSYDGSGRGFFDFWGQGTLVTVSS"
XX
SQ   Sequence 840 BP; 157 A; 227 C; 265 G; 191 T; 0 other;
     atgtggtggc gcctgtggtg gctgctgctg ctgctgctgc tgctgtggcc catggtgtgg        60
     gccgatattg tgctgacgca gtctccaggc accctgtctt tgtctgcagg ggaaagagcc       120
     accctctcct gcagggccag tcagagtgtt agcagcggct ccttagcctg gtaccagcag       180
     aaacctggtc aggctcccag gctcctcatc tacggtgcat ccaccagggc cactggcatc       240
     ccagacaggt tcagtggcag tgggtctggg acagacttca ctctcacaat cggcagactg       300
     gagcctgaag atctcgcagt atattactgt cagcagtatg gtacctcacc gtacactttt       360
     ggccagggga ccaaagtgga tatcaaacgt ggtggcggtg gctcgggcgg tggcggttca       420
     ggtggcggtg gctctagatc ttcccaggtc cagcttgtgc agtctggggc tgaggtgaag       480
     aagcctgggt cctcggtgca ggtctcctgc aaggcctctg gaggcacctt cagcatgtat       540
     ggtttcaact gggtgcgaca ggcccctgga catggccttg agtggatggg agggatcatc       600
     cctatctttg gtacatcaaa ctacgcacag aagttccggg gcagagtcac gtttaccgcg       660
     gaccaagcca cgagcacagc ctacatggag ctgaccaacc tgcgatctga cgacacggcc       720
     gtctattatt gtgcgagaga ttttggcccc gactgggaag acggtgattc ctatgatggt       780
     agtggccggg ggttctttga cttctggggc cagggaaccc tggtcaccgt ctcctcatga       840
//
"""
    expected_output = ['AC   JA013192;', 'FT                   /organism="Homo sapiens"', 
                       'FT                   /protein_id="CCA61653.1"', 
                       'FT                   /translation="MWWRLWWLLLLLLLLWPMVWADIVLTQSPGTLSLSAGERATLSCR',
                       'FT                   ASQSVSSGSLAWYQQKPGQAPRLLIYGASTRATGIPDRFSGSGSGTDFTLTIGRLEPED', 
                       'FT                   LAVYYCQQYGTSPYTFGQGTKVDIKRGGGGSGGGGSGGGGSRSSQVQLVQSGAEVKKPG', 
                       'FT                   SSVQVSCKASGGTFSMYGFNWVRQAPGHGLEWMGGIIPIFGTSNYAQKFRGRVTFTADQ', 
                       'FT                   ATSTAYMELTNLRSDDTAVYYCARDFGPDWEDGDSYDGSGRGFFDFWGQGTLVTVSS"']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compile.....')
    parser.add_argument(
        '--emblfile', help='EMBL file', required=True)
    parser.add_argument(
        '--newfastadir', help='Directory for location of new fasta files', required=True)
    args = parser.parse_args()

    get_data4fasta(args.emblfile, args.newfastadir)