#!/usr/bin/env python3

# receive file in Illumina SNP_Map format, and convert it into MAP file (for PLINK)
###################################################################################
# Written by Moran Neuhof, 2017

from sys import argv, exit

def fix_sex_chrs(chr_name):
    """Fix the chr names to fit the PLink standard. 
    Return fixed chrname"""
    chr_name_dict = {'X': '19',
                     'Y': '20'}  # must replace sex chromosomes with numbers, for Plink

    try:
        return chr_name_dict[chr_name]  # if it's a sex chromosome
    except KeyError:
        return chr_name  # if it is an autosome



if len(argv) < 2:
    print("Error: Missing map.txt file!")
    exit(1)
else:
    SNP_map_file = argv[1]  # use input

# set output file name:
SNP_map_output = SNP_map_file + ".ped"


with open(SNP_map_file, 'r') as infile:
    split_lines = (line.split('\t') for line in infile)  # splitting each line
    write_line_list = [[fix_sex_chrs(split_line[2]), split_line[1], '0', split_line[3]] for split_line in split_lines][1:]
    # removed the first line

# writing the results to a map file:
with open(SNP_map_output, 'w') as outfile:
    for line_to_write in write_line_list:
        print("\t".join(line_to_write), end="\n", file=outfile)

print(f"Done, output map file saved to {SNP_map_output}")
