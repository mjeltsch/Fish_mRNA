#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This program takes one or multiple nucleotide sequence fasta files
# (manually exported from the PhyloFish database, http://phylofish.sigenae.org/)
# and extracts the ORFs longer than "minimal_length". Then, it
# checks for the presence of the PDGF/VEGF homology domain,
# using the RE pattern "pattern". It extracts all the amino acid
# sequence of the longest reading frame

# This program has been tested on Ubuntu 18.04 and 20.04.
# It has a bunch of requirements, the most important being:
# 
# Biopython, t_coffee, mpirun, phyml, python modules argparse, re, shutil
#
# Usage: ./extract_ORF.py fasta_mRNA_contigs/*.fasta

import argparse,re,shutil,os
from Bio import SeqIO
from phylolib import execute_subprocess

def run():
    # Relaxed PDGF signature
    pattern = 'P.?C.{2,8}C.?G.?C' # 17 from 61
    #pattern = 'P.?C.{3,7}C.?G.?C' # 16 from 61
    #pattern = 'P.C.{4,6}C.G.C' #16 from 61
    minimal_length = 300
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfiles', nargs='+')
    args = parser.parse_args()
    i = 0
    k = 0
    for inputfile in args.inputfiles:
        base_name = os.path.basename(inputfile)
        aa_sequence_list = []
        j = 0
        records = SeqIO.parse(inputfile, 'fasta')
        for record in records:
            i += 1
            j += 1
            # Check both fwd and reverse strands
            for strand, seq in (1, record.seq), (-1, record.seq.reverse_complement()):
                # Check all three reading frames
                for frame in range(3):
                    # Division // returns ganzzahligen Anteil
                    length = 3 * ((len(seq)-frame) // 3)
                    # Take all peptides except the last
                    for pro in seq[frame:frame+length].translate(table = 1).split("*")[:-1]:
                        # Look for start codon (ATG/M)
                        if 'M' in pro:
                            orf = pro[pro.find('M'):]
                            pos = seq[frame:frame+length].translate(table=1).find(orf)*3 + frame +1
                            if len(orf)*3 +3 > minimal_length:
                                # Search for PDGF/VEGF pattern (for the first match ONLY)
                                if re.search(pattern, str(orf)):
                                    print("{}...{} - length {}aa, strand {}, frame {}, pos {}, name {}".format(orf[:3], orf[-3:], len(orf)+1, strand, frame, pos, record.id))
                                    aa_sequence_list.append(orf)
                                    k += 1
                                else:
                                    print('No ORFs found > {0} aa, matching {1} found.'.format(minimal_length, pattern))
    
        # Write all ORFs into one fasta file
        outputfile = 'results/' + base_name + '_aa.fa'
        with open(outputfile, "w") as output_handle:
            for sequence in aa_sequence_list: 
                output_handle.write('>' + record.id + '\n' + str(sequence) + '\n')
        print('PDGF signature found in {0} of {1} sequences from species {2}.'.format(len(aa_sequence_list), str(j), record.id.split('_')[0]))
        # Add Danio rerio sequences to fasta files
        mergefile = 'results/' + base_name + '_merged.fa'
        with open(mergefile,'wb') as wfd:
            for file in [outputfile, 'All_Danio_rerio.fa']:
                with open(file,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)

        # Do alignment
        aligned_file = base_name + '_aligned.fa'
        mergefile = base_name + '_merged.fa'
        comment = 'Aligning...'
        command = 't_coffee -mode mcoffee -in {0} > {1}'.format(mergefile, aligned_file)
        execute_subprocess(comment, command, working_directory = 'results')

        PROCESSORS = 6
        # Use -1 for testing, and something between 100-1000 for production use
        BT_REPLICATES = -1
        aligned_phy_file = base_name + '_aligned.phy'
        codename_file = base_name + '_codenames.txt'
        aligned_file_encoded =  base_name + '_encoded.fa'
        treefile_encoded = base_name + '_aligned.phy_phyml_tree.txt'
        treefile_decoded = base_name + '_decoded.phy_phyml_tree.txt'

        # Encode sequence names because phyml butchers them
        comment = 'Generating codename file...'
        command = 't_coffee -other_pg seq_reformat -in {0} -output code_name > {1}'.format(aligned_file, codename_file)
        execute_subprocess(comment, command, working_directory = 'results')

        comment = 'Encoding sequence names...'
        command = 't_coffee -other_pg seq_reformat -code {0} -in {1} > {2}'.format(codename_file, aligned_file, aligned_file_encoded)
        execute_subprocess(comment, command, working_directory = 'results')

        # Convert fasta alignment into phylip alignment (beecause that what phyml understands)
        comment = 'Converting alignment from fasta to phylip format...'
        command = 't_coffee -other_pg seq_reformat -in {0} -output phylip > {1}'.format(aligned_file_encoded, aligned_phy_file)
        execute_subprocess(comment, command, working_directory = 'results')

        # FINALLY: Build the tree
        comment = 'Treebuilding...'
        command = 'mpirun -n {0} phyml-mpi --no_memory_check -s BEST -i {1} -d aa -b {2}'.format(PROCESSORS, aligned_phy_file, BT_REPLICATES)
        execute_subprocess(comment, command, working_directory = 'results')

        # Decode sequence names in the main tree file that humans can interpret the result 
        comment = 'Decoding sequence names...'
        command = 't_coffee -other_pg seq_reformat -decode {0} -in {1} > {2}'.format(codename_file, treefile_encoded, treefile_decoded)
        execute_subprocess(comment, command, working_directory = 'results')

    print('PDGF signature found in {0} of {1} sequences.'.format(str(k), str(i)))

if __name__ == '__main__':
    run()
