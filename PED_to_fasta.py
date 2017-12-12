#!/usr/bin/env python3

# receive a PED file (PLINK textual format), and convert it into a fasta file
###############################################################################
# Written by Moran Neuhof, 2017

import os
import random
from sys import argv,exit


def compare_nuc_couple(nuc1, nuc2, ambiguous):
    """Receive two nucleotides and returns one consensus nucleotide."""

    ambiguous_conversion_dict = {('A', 'C'): 'M',
                                 ('C', 'A'): 'M',
                                 ('A', 'G'): 'R',
                                 ('G', 'A'): 'R',
                                 ('A', 'T'): 'W',
                                 ('T', 'A'): 'W',
                                 ('G', 'C'): 'S',
                                 ('C', 'G'): 'S',
                                 ('C', 'T'): 'Y',
                                 ('T', 'C'): 'Y',
                                 ('G', 'T'): 'K',
                                 ('T', 'G'): 'K'}

    if nuc1 == nuc2:        # if they are the same
        if nuc1 == '0':     #
            cons_nuc = 'N'  # uncovered
        else:
            cons_nuc = nuc1
    else:  # ambiguous
        if ambiguous:
            cons_nuc = ambiguous_conversion_dict[(nuc1, nuc2)]
        else: 
            cons_nuc = 'N'
    return cons_nuc


def compare_nuc_couple_no_ambig(nuc1, nuc2):
    """Receive two nucleotides and returns one consensus nucleotide."""

    if nuc1 == nuc2:
        if nuc1 == '0':  # if they are the same
            return 'N'  # uncovered
        else:
            return nuc1
    else:  # ambiguous
        return 'N'


def parse_ped_line(ped_line, with_ambiguous):
    """Compare every two nucleotides in each line."""
    
    # each 2 nucleotides are actually one SNP
    ped_line_pos1 = [x for i, x in enumerate(ped_line) if i % 2 == 0]
    ped_line_pos2 = [x for i, x in enumerate(ped_line) if i % 2 == 1]

    # iterate over each homologue nucleotides in the same position
    if with_ambiguous:
        nuc_list = list(map(compare_nuc_couple, ped_line_pos1, ped_line_pos2))  # ambiguous chars
    else:
        nuc_list = list(map(compare_nuc_couple_no_ambig, ped_line_pos1, ped_line_pos2))  # just AGCT or N
    return "".join(nuc_list)



def choose_outgroup(outgroup_list):
    """Returns a random representative of the outgroups"""
    return random.choice(outgroup_list)


def create_outfile_name(infile_name, ambig):
    """Receive file name and ambig (Boolean) and returns an output file name (string)"""
    # add ambig to filename
    if ambig == True:
        ambig = '.ambig'
    else:
        ambig = ''

    if infile_name.endswith('.ped'):
        fasta_outfile = infile_name.replace('.ped', f'{ambig}.fasta')
    elif infile_name.endswith('.PED'):
        fasta_outfile = infile_name.replace('.PED', f'{ambig}.fasta')
    else:
        fasta_outfile = infile_name + f'_{ambig}.fasta'
    return fasta_outfile


# check input
if len(argv) < 2:
    print("Error: Missing PED file!")
    print("Usage:")
    print(f"{argv[0]} <PED file> (--ambig)")
    print(f"--ambig: The flag includes ambiguous nucleotides. If not present, ambiguous positions will be turned into N's.")
    exit(1)
else:
    ped_file = argv[1]  # use input
    if len(argv) > 2 and argv[2] in ('--ambig', '--ambiguous', '-ambig', '-ambiguous'):
        AMBIG = True  # set ambiguous mode to includ ambiguous
        print("Using ambiguous characters.")
    else:
        AMBIG = False
        print("Not using ambiguous characters.")

fasta_outfile = create_outfile_name(ped_file, AMBIG)  # output

# The following part sets specific outgroups for the data. 
# You probably won't need it 
#########################################################
# # outgroups
# BABA_outgroups = [f'>BABA_BABA_{i+1}_W' for i in range(4)]
# PHAF_outgroups = [f'>PHAF_PHAF_{i+1}_W' for i in range(8)]
# SBSB_outgroups = [f'>SBSB_SBSB_{i+1}_W' for i in range(11)]
# SCEL_outgroups = [f'>SCEL_SCEL_{i+1}_W' for i in range(6)]
# SVSV_outgroups = [f'>SVSV_SVSV_{i+1}_W' for i in range(10)]
# # a list of all outgroups
# all_outgroups = [*BABA_outgroups, *PHAF_outgroups, *SBSB_outgroups, *SCEL_outgroups, *SVSV_outgroups]
# # a list of the chosen outgroups to keep
# outgroups_to_keep = [choose_outgroup(outgroup_list) for outgroup_list in
#                      (BABA_outgroups, PHAF_outgroups, SBSB_outgroups, SCEL_outgroups, SVSV_outgroups)]
#

with open(ped_file, 'r') as infile:
    with open(fasta_outfile, 'w') as outfile:
        for line in infile:
            split_line = line.strip().split(' ')  # split every line
            header = f">{split_line[0]}_{split_line[1]}"  # the fasta header
            
            # again, this part is about outgroups, so I commented it out:
            
            # check if it's outgroup and we need to keep it
            # if header in all_outgroups and header not in outgroups_to_keep:
            #     print(f"Skipping outgroup {header}")
            #     continue # if it is an outgroup that we didn't choose

            ped_sequence = split_line[6:]
            fasta_line = parse_ped_line(ped_sequence, with_ambiguous=AMBIG)
            print(f"{header}\n{fasta_line}", end='\n', file=outfile)  # writing to fasta file

# print(f"Chosen outgroups: {','.join(x[1:] for x in outgroups_to_keep)}")
print(f"Done, output fasta file saved to {fasta_outfile}")
