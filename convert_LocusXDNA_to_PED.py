#!/usr/bin/env python3

# receive file in Illumina LocusXDNA format, and convert it into PED
####################################################################
# Written by Moran Neuhof, 2017

from sys import argv, exit


def parse_XDNA_file(XDNA_file):
    """Receive a XDNA filename,
    return a locus_list (basically a numbered list),
           a SNP_list with the same length
           a dictionary with a list of the SNP values in each sample"""

    genotype_dict = {}  # initialize

    with open(XDNA_file, 'r') as infile:
        for split_line in (line.strip().split(',') for line in infile):
            line_len = len(split_line)
            if line_len > 2 and split_line[1] == 'snps':
                SNP_list = split_line[3:]
            elif line_len > 2 and split_line[1] == 'locus':
                locus_list = split_line[3:]
            elif line_len > 6 and split_line[0] != 'PorcineSNP60v2_B' and split_line[6] == 'calls':  # then this is a genotyping line
                genotype_dict[split_line[0]] = split_line[8:]  # saving the genotype for each animal
    return locus_list, SNP_list, genotype_dict


def convert_real_val_to_two_explicit_bases(A, B, real_val):
    """Receives the value of A, B for the position and the real_val for the sample in that position, and returns it as
    a string for writing in ped"""

    convert_dict = {'A': [A, A],
                    'B': [B, B],
                    'H': [A, B],
                    'U': ['0', '0']}  # fixing the phased data
    return " ".join(convert_dict[real_val])



# convert to ped format
def convert_genotype_to_ped(genotype_list, SNP_value_list):

    split_genotype_ref_values = (x.split('/') for x in genotype_list)
    return [convert_real_val_to_two_explicit_bases(A, B, real_val) for (A, B), real_val in
                        zip(split_genotype_ref_values, SNP_value_list)]


def fix_sample_name(sample_name):

    sample_name_dict = {'MM10': 'MM010',
                        'MM11': 'MM011',
                        'MM14': 'MM014'}  # had to fix the sample names

    try:
        return sample_name_dict[sample_name]  # if it's a sex chromosome
    except KeyError:
        return sample_name  # if it is an autosome


def find_sex(sample):
    """Receives a sample name and return its sex (according to PLINK documentation):
    1 = Male
    2 = Female
    0 = Unknown"""
    sex_dict = {'MM010': '2',
                'MM011': '1',
                'MM280': '2',
                'MM281': '1',
                'MM282': '2'}
    # all data in the dictionary taken from the file "Modern pigs that were sent to SNPs.xlsx"
    # by Meirav Meiri

    try:
        return sex_dict[sample]  # return its sex
    except KeyError:  # no sex for this poor pig
        return '0'


# Parse the file
# Locus_XDNA_filename
if len(argv) < 2:
	print("Error: Missing Locus_XDNA file!")
	exit(1)
else:
	Locus_XDNA_filename = argv[1]  # use input

# set output file name:
out_filename = Locus_XDNA_filename + ".ped"


locus_list, SNP_list, genotype_dict = parse_XDNA_file(Locus_XDNA_filename)

with open(out_filename, 'w') as outfile:
    for sample_name, SNP_values_per_sample in genotype_dict.items():
        ped_SNP_val_list = convert_genotype_to_ped(SNP_list, SNP_values_per_sample)  # get a list of PED values
        print(len(ped_SNP_val_list))
        sample_name = fix_sample_name(sample_name)  # fixed the three problematic samples
        # create a list of the first 6 ped values:
        # Family ID, IID, father, mother, sex, phenotype
        sample_header_list = ['Israel', sample_name, "0", "0", find_sex(sample_name), "0"]  # change to fit other samples
        print(" ".join(sample_header_list))
        ped_line = f"{' '.join(sample_header_list)} {' '.join(ped_SNP_val_list)}"
        print(len(ped_line))
        print(ped_line, end='\n', file=outfile)

print(f"Done.\nPED file created in {out_filename}")