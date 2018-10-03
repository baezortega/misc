#!/usr/bin/env python
#coding: utf-8

# TransductionFinder.py
# Adrian Baez-Ortega
# 2018


"""
Given the paths to a set of sample and matched control BAM files (in a tab-delimited 
file with two columns), and a set of target region coordinates (in BED format), this 
script identifies retrotransduction events in each sample, by searching for discordant 
read pairs (in which one of the reads is mapped to a distant position or a different 
chromosome) in each of the input regions. The script produces a tab-delimited file 
containing, for each candidate event:
  the name of the sample and matched control BAM files;
  the transduction type (interchromosomal or intrachromosomal);
  the coordinates of the input source region; 
  the coordinates of the region covered by reads supporting the event at the source locus; 
  the coordinates of the region covered by reads supporting the event at the insertion locus; 
  the number of read pairs supporting the event in the sample;
  the number of read pairs supporting the event in the matched control.
(NB. This script requires Python 2.7 and the Python library Pysam version >=0.9.1.4.)
"""


import sys, os, argparse, time, pysam


####################
# Argument parsing #
####################
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='Identifies retrotransduction events in a set of genomic regions in BAM files')
parser.add_argument('sample_list', help='Tab-delimited text file with two columns, containing the paths to sample and matched control BAM files, one file pair per line')
parser.add_argument('regions_file', help='BED file containing coordinates of the target source genomic regions')
parser.add_argument('out_file', help='Path to output file (will be overwritten if exists)')
parser.add_argument('-r', '--reads', type=int, default=3, help='Minimum number of read pairs supporting a transduction in the sample BAM')
parser.add_argument('-d', '--dist', type=int, default=100, help='Minimum distance (in kb) for intra-chromosomal transduction')

args = parser.parse_args()
SAMPLE_LIST = args.sample_list
REGIONS_FILE = args.regions_file
OUT_FILE = args.out_file
MIN_READS = args.reads
MIN_DIST = args.dist
MIN_DIST_BP = MIN_DIST * 1000

# List of valid chromosome names
VALID_CHROMS = [str(x) for x in range(1, 23)] + ['X', 'Y']


###############
# Main script #
###############

print '\nRunning TransductionFinder.py on ' + time.asctime(time.localtime(time.time()))
print '\nWith arguments:'
print 'Sample and control BAM file list:     ', SAMPLE_LIST
print 'Source region coordinates file:       ', REGIONS_FILE
print 'Output file:                          ', OUT_FILE
print 'Min. supporting read pairs in sample: ', MIN_READS
print 'Min. intrachromosomal distance:       ', MIN_DIST, 'kb'


# Process each sample
with open(SAMPLE_LIST, 'r') as samples, open(OUT_FILE, 'w') as out:

    # Write output file header
    out.write('## Output by TransductionFinder.py on ' + time.asctime(time.localtime(time.time())) + '\n' + \
              '## Input BAM file list: ' + SAMPLE_LIST + '\n' + \
              '## Input BED file: ' + REGIONS_FILE + '\n' + \
              '## Minimum supporting read pairs per event in sample: ' + str(MIN_READS) + '\n' + \
              '## Minimum distance between intrachromosomal read clusters: ' + str(MIN_DIST) + ' kb\n' + \
              '## Output fields:\n' + \
              '##  sample_name: name of sample BAM file (without extension)\n'
              '##  control_name: name of matched control BAM file (without extension)\n'
              '##  trands_type: transduction type (interchromosomal or intrachromosomal)\n' + \
              '##  input_source_chrom: chromosome of input source region\n' + \
              '##  input_source_start: start position of input source region\n' + \
              '##  input_source_end: end position of input source region\n' + \
              '##  supp_source_start: start position of supported region at source\n' + \
              '##  supp_source_end: end position of supported region at source\n' + \
              '##  supp_dest_chrom: chromosome of supported region at destination\n' + \
              '##  supp_dest_start: start position of supported region at destination\n' + \
              '##  supp_dest_end: end position of supported region at destination\n' + \
              '##  supp_reads_sample: number of read pairs supporting transduction in the sample BAM\n' + \
              '##  supp_reads_control: number of read pairs supporting transduction in the matched control BAM\n' + \
              '## (supported regions are defined by the START positions of supporting reads)\n' +\
              '\t'.join(['sample_name', 'control_name',
                         'transd_type', 
                         'input_source_chrom', 'input_source_start', 'input_source_end',
                         'supp_source_start', 'supp_source_end',
                         'supp_dest_chrom', 'supp_dest_start', 'supp_dest_end', 
                         'supp_reads_sample', 'supp_reads_control']) + '\n')
    
    for line in samples:
        sample_path, control_path = line.strip().split('\t')
        sample_name = os.path.splitext(os.path.basename(sample_path))[0]
        control_name = os.path.splitext(os.path.basename(control_path))[0]
        
        print '\nProcessing sample', sample_name
        
        # Open BAM files
        sample_bam = pysam.AlignmentFile(sample_path, "rb")
        control_bam = pysam.AlignmentFile(control_path, "rb")
        
        # Process input regions
        with open(REGIONS_FILE, 'r') as regions:            
            for region in regions:
                [chrom, start, end] = region.strip().split('\t')
                region_coords = chrom + ':' + start + '-' + end
                
                # Initialise candidate supporting read pair sets
                # Candidate read pairs are clustered in groups with max inter-mate distance MIN_DIST_BP
                # Each candidate read pair group follows structure: [[read_pos1, ..., read_posN], 
                #                                                    [mate_pos1, ..., mate_posN]]
                cand_inter = {}
                cand_intra = []
                
                print '  Extracting reads in', region_coords
                
                # Retrieve alignments
                reads = sample_bam.fetch(chrom, int(start), int(end))
                for read in reads:
                
                    # Use only mapped, paired reads which are not PCR duplicates
                    if not (read.is_duplicate or read.is_unmapped or read.mate_is_unmapped):
                        
                        # Retrieve read and mate coordinates
                        read_pos = read.reference_start
                        mate_chrom = read.next_reference_name
                        mate_pos = read.next_reference_start
                        
                        # Identify interchromosomal candidates
                        if mate_chrom != chrom:
                            if mate_chrom in VALID_CHROMS:
                                if mate_chrom in cand_inter:
                                    # Check whether mate clusters with previously detected mates
                                    found = False
                                    for cand in cand_inter[mate_chrom]:
                                        # If distance between mate_pos and mates in group < MIN_DIST_BP
                                        if abs(mate_pos - min(cand[1])) < MIN_DIST_BP or \
                                           abs(mate_pos - max(cand[1])) < MIN_DIST_BP:
                                            # Add read pair to group
                                            cand[0].append(read_pos)
                                            cand[1].append(mate_pos)
                                            found = True
                                            break
                                    # Otherwise, create new group
                                    if not found:
                                        cand_inter[mate_chrom].append([[read_pos], [mate_pos]])
                            
                                else:
                                    cand_inter[mate_chrom] = [ [[read_pos], [mate_pos]] ]
                        
                        # Identify intrachromosomal candidates
                        elif abs(mate_pos - read_pos) > MIN_DIST_BP:
                            # Check whether mate clusters with previously detected mates
                            found = False
                            for cand in cand_intra:
                                # If distance between mate_pos and mates in group < MIN_DIST_BP
                                if abs(mate_pos - min(cand[1])) < MIN_DIST_BP or \
                                   abs(mate_pos - max(cand[1])) < MIN_DIST_BP:
                                    # Add read pair to group
                                    cand[0].append(read_pos)
                                    cand[1].append(mate_pos)
                                    found = True
                                    break
                            # Otherwise, create new group
                            if not found:
                                cand_intra.append([[read_pos], [mate_pos]])
                
                
                # Output candidate read groups if there is enough support
                # Intrachromosomal
                type = 'intrachrom'
                for cand in cand_intra:
                    if len(cand[0]) >= MIN_READS:
                    
                        source_start = min(cand[0])
                        source_end = max(cand[0])
                        dest_chrom = chrom
                        dest_start = min(cand[1])
                        dest_end = max(cand[1])
                        supp_pairs_sample = len(cand[0])
                        
                        # Count supporting read pairs in matched control
                        supp_pairs_control = 0
                        reads = control_bam.fetch(chrom, source_start - 200, source_end + 200)
                        for read in reads:
                            if not (read.is_duplicate or read.is_unmapped or read.mate_is_unmapped):                        
                                mate_chrom = read.next_reference_name
                                mate_pos = read.next_reference_start
                                
                                if mate_chrom == dest_chrom and \
                                   mate_pos >= dest_start - 200 and mate_pos <= dest_end + 200:
                                    supp_pairs_control += 1
                        
                        # Output
                        out.write('\t'.join([sample_name, control_name,
                                             type, chrom, start, end,
                                             str(source_start), str(source_end),
                                             dest_chrom, str(dest_start), str(dest_end),
                                             str(supp_pairs_sample), str(supp_pairs_control)]) + '\n')
                        
                        source_coords = chrom + ':' + str(source_start) + '-' + str(source_end)
                        dest_coords = dest_chrom + ':' + str(dest_start) + '-' + str(dest_end)
                        print '    Intrachromosomal candidate: ' + source_coords + ' to ' + dest_coords + \
                              ' (' + str(supp_pairs_sample) + ' read pairs in sample, ' + str(supp_pairs_control) + ' read pairs in control)'
                
                # Interchromosomal
                type = 'interchrom'
                for dest_chrom, cands in sorted(cand_inter.items()):
                    for cand in cands:
                        if len(cand[0]) >= MIN_READS:
                        
                            source_start = min(cand[0])
                            source_end = max(cand[0])
                            dest_start = min(cand[1])
                            dest_end = max(cand[1])
                            supp_pairs_sample = len(cand[0])
                        
                            # Count supporting read pairs in matched control
                            supp_pairs_control = 0
                            reads = control_bam.fetch(chrom, source_start - 200, source_end + 200)
                            for read in reads:
                                if not (read.is_duplicate or read.is_unmapped or read.mate_is_unmapped):                        
                                    mate_chrom = read.next_reference_name
                                    mate_pos = read.next_reference_start
                                
                                    if mate_chrom == dest_chrom and \
                                       mate_pos >= dest_start - 200 and mate_pos <= dest_end + 200:
                                        supp_pairs_control += 1
                            
                            # Output
                            out.write('\t'.join([sample_name, control_name,
                                                 type, chrom, start, end,
                                                 str(source_start), str(source_end),
                                                 dest_chrom, str(dest_start), str(dest_end),
                                                 str(supp_pairs_sample), str(supp_pairs_control)]) + '\n')
                            
                            source_coords = chrom + ':' + str(source_start) + '-' + str(source_end)
                            dest_coords = dest_chrom + ':' + str(dest_start) + '-' + str(dest_end)
                            print '    Interchromosomal candidate: ' + source_coords + ' to ' + dest_coords + \
                                  ' (' + str(supp_pairs_sample) + ' read pairs in sample, ' + str(supp_pairs_control) + ' read pairs in control)'

        sample_bam.close()
        control_bam.close()

print '\nDone\n'
