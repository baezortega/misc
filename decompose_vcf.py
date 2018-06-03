#!/usr/bin/env python

# decompose_vcf.py
# Adrian Baez-Ortega, 2018

# This script has been tested in Python 2.7 in Unix systems


"""
This script processes a VCF, or compressed (gzipped) VCF file, extracting
the variant metadata fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO),
and those values specified by the user (by default, all) from the genotype
information of each variant in each sample. Each of these values is output
as a table to a separate output text file.
"""


import sys
import gzip
import os
import re


# Constant: position of FORMAT field in VCF
FORMAT_COL = 8


# Helper function: open a VCF or gzipped VCF file
def open_vcf(vcf_path, gzipped):
    if gzipped:
        return open(vcf_path, 'r')
    else:
        return gzip.open(vcf_path, 'r')


# If not 1 argument: print help
if len(sys.argv) < 2:
    print '\ndecompose_vcf.py'
    print 'This script processes a VCF, or compressed (gzipped) VCF file, extracting'
    print 'the variant metadata fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO),'
    print 'and those values specified by the user (by default, all) from the genotype'
    print 'information of each variant in each sample. Each of these values is output'
    print 'as a table to a separate output text file.'
    print '\nInput: Path to input VCF (.vcf) or gzipped VCF (.vcf.gz) file.'
    print '       [Optional] List of value names of interest (as shown in the FORMAT field)'
    print '\nUsage: decompose_vcf.py /path/to/file.vcf.gz [GT GL NR NV ...]\n'
    sys.exit(0)


script, vcf_path = sys.argv[0:2]
if len(sys.argv) > 2:
    value_names = sys.argv[2:]


# Check if VCF is gzipped
if vcf_path[-7:] == '.vcf.gz':
    ext_len = 7
    gzipped = True
    
elif vcf_path[-4:] == '.vcf':
    ext_len = 4
    gzipped = False
    
else:
    print '\nERROR: Invalid file extension. Only extensions .vcf and .vcf.gz are accepted.\n'
    sys.exit(1)


# If value names are specified, check that they appear in the FORMAT field
# If no value names are specified, retrieve all names from FORMAT field
with open_vcf(vcf_path, 'r') as vcf:
    line = vcf.readline()
    while line.startswith('#'):
        line = vcf.readline()
    
    format_names = line.split('\t')[FORMAT_COL].split(':')
    if len(sys.argv) > 2:
        for name in value_names:
            if name not in format_names:
                print '\nERROR: Name', name, 'not found in FORMAT field of the VCF file.\n'
                sys.exit(1)
    else:
        value_names = line.split('\t')[FORMAT_COL].split(':')

    # Match specified value names to their positions in FORMAT
    value_idx = []
    for name in value_names:
        value_idx.append([i for i in range(len(format_names)) if format_names[i] == name][0])


# Compose output file paths
out_files = []
for name in value_names:
    out_files.append(vcf_path[:-ext_len] + '_' + name + '.txt')
out_files.append(vcf_path[:-ext_len] + '_Meta.txt')


# Print info
print '\nInput file:'
print '  ' + vcf_path
print 'Extracting variant metadata and values:' 
print '  ' + ', '.join(value_names)
print 'Output files:'
for of in out_files:
    print '  ' + of


# Process VCF
print 'Processing VCF...'
try:

    # Open all output files as a list of file descriptors
    out = [open(of, 'w') for of in out_files]

    with open_vcf(vcf_path, 'r') as vcf:
        for line in vcf:
            if not line.startswith('##'):
                
                cols = line.strip().split('\t')
                
                # Copy column header to output files
                if 'POS' in cols[1]:
                    header_meta = cols[:FORMAT_COL]
                    header_samples = cols[FORMAT_COL+1:]
                    out[-1].write('\t'.join(header_meta) + '\n')
                    for of in out[:-1]:
                        of.write('\t'.join(header_samples) + '\n')
                
                else:
                    # Output metadata
                    metadata = cols[:FORMAT_COL]
                    out[-1].write('\t'.join(metadata) + '\n')
                    
                    # Extract selected values from sample genotypes
                    values = [[] for name in value_names]
                    data = cols[FORMAT_COL+1:]
                    for field in data:
                        fld = field.split(':')
                        for i in range(len(value_idx)):
                            values[i].append(fld[value_idx[i]])
                    
                    # Output extracted values
                    for i in range(len(values)):
                        out[i].write('\t'.join(values[i]) + '\n')

finally:
    for of in out:
        of.close()

print 'Done!\n'
