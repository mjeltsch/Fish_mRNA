# Fish_mRNA
This script looks through mRNA sequences and extracts PDGF/VEGF-like ORFs 

This program takes one or multiple nucleotide sequence fasta files (manually exported from the PhyloFish database, http://phylofish.sigenae.org/) and extracts the ORFs longer than "minimal_length". Then, it checks for the presence of the PDGF/VEGF homology domain, using the RE pattern "pattern". It extracts all the amino acid sequence of the longest reading frame.

This program has been tested on Ubuntu 18.04 and 20.04.
It has a bunch of requirements, the most important being:
 
- Biopython
- t_coffee
- mpirun
- phyml
- python modules argparse, re, shutil

Usage:

    ./extract_ORF.py fasta_mRNA_contigs/*.fasta
