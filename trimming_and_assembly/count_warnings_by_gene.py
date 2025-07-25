#!/usr/bin/python

import sys
import os.path

gene_list = []
sample_list = []

def read_gene_file():
    gene_file = sys.argv[1]
    with open(gene_file) as gf_obj:
        for line in (gf_obj):
            gene_list.append(line.rstrip())
    #print(gene_list[0:9])
    
def read_sample_file():
    sample_file = sys.argv[2]
    with open(sample_file) as sf_obj:
        for line in (sf_obj):
            sample_list.append(line.rstrip())


def read_warnings_by_sample(sample_code,warning):
    file_suffixes = {
    'chimeric_stitched' : '_genes_derived_from_putative_chimeric_stitched_contig.csv',
    'internal_stop' : '_genes_with_non_terminal_stop_codons.txt',
    'stitched_contig' : '_genes_with_stitched_contig.csv',
    'failed_spades' : 'failed_spades.txt'
    }
    if warning == "failed_spades":
        warning_file=f"{sample_code}/{file_suffixes[warning]}"
    else:
        warning_file = f"{sample_code}/{sample_code}{file_suffixes[warning]}"
    #print(warning_file)
    gene_warnings = {}
    if os.path.isfile(warning_file):
        with open(warning_file) as wf_obj:
            for line in (wf_obj):
                gene_warnings[line.rstrip()] = 1
    return gene_warnings


read_gene_file()
read_sample_file()

print("sample",end="")
for gene in gene_list:
    print("\t",gene,end="")
print()

for sample_code in sample_list:
    warnings=read_warnings_by_sample(sample_code,'stitched_contig')
    print(sample_code,end="")
    for gene in gene_list:
        if gene in warnings.keys():
            print("\t","1",end="")
        else:
            print("\t","0",end="")
    print()
    
