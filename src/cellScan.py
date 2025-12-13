#!/usr/bin/env python3

# This script performs low-level scanning of single-cell BAM reads
# to evaluate whether reads support reference or alternative alleles
# at known SNV positions. It operates at the read–sequence level,
# explicitly matching sequence motifs around SNV sites.

import argparse
import sys
import os
import logging
import shutil
import glob
import re
import pysam
import time
import subprocess
import pandas as pd
import numpy as np
import gzip
from pysam import VariantFile
import multiprocessing as mp
from multiprocessing import Pool



def robust_get_tag(read, tag_name):
    # Safely retrieve a BAM tag from a read.
    # If the tag does not exist, return "NotFound" instead of raising an error.
    # This is useful for optional tags such as cell barcodes (CB).
    try:
        return read.get_tag(tag_name)
    except KeyError:
        return "NotFound"



def rev_compl(st):
    # Compute the reverse complement of a DNA sequence.
    # Used when evaluating motifs on the opposite strand.
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))


# Hard-coded reference inputs (likely used for testing / debugging)
# Reference genome FASTA
in_fasta = "/rsrch3/scratch/bcb/jdou1/scAncestry/ref/fasta/genome.fa"

# Input VCF containing SNV sites
in_vcf = "/rsrch3/scratch/bcb/jdou1/bam_process/chr1.gp.vcf"

# Input BAM containing aligned single-cell reads
in_bam = "/rsrch3/scratch/bcb/jdou1/bam_process/chr1.filter.targeted.bam"


# Open reference genome and VCF
ref_fa = pysam.FastaFile(in_fasta)
vcf = VariantFile(in_vcf)

# Dictionary storing SNV metadata indexed by an integer counter
snv_info = {}
index = 0

# Length of motif on each side of the SNV
# Total motif length = 2 * motif_len + 1
motif_len = 3


# Iterate through SNVs in the VCF and extract sequence context
for rec in vcf.fetch():
    # Extract reference sequence surrounding the SNV
    # (3 bp upstream + SNV + 3 bp downstream)
    seq = ref_fa.fetch(rec.chrom, rec.pos - 4, rec.pos + 3)

    # Reverse complement of the motif (not directly used later)
    seq_rev_compl = rev_compl(seq)

    # Unique SNV identifier
    id = str(rec.chrom) + ":" + str(rec.pos) + ":" + rec.ref + ":" + rec.alts[0]

    # Store SNV information
    mydict = dict(
        chr=rec.chrom,
        pos=rec.pos,
        motif_pos=seq,
        ref_allele=rec.ref,
        alt_allele=rec.alts[0]
    )

    index += 1
    snv_info[index] = mydict


# Total number of SNVs
snv_tol = index

# Initialize scanning indices
index = 1
lower_index = 1


# Counters for read statistics
read_tol = 0          # total reads processed
read_cover = 0        # reads covering at least one SNV
read_wild = 0         # reads supporting reference allele
read_mutated = 0      # reads supporting alternative allele
read_noAllele = 0     # reads covering SNV but matching neither allele

# Open BAM file
infile = pysam.AlignmentFile(in_bam, "rb")

# Log file for ambiguous cases
fp = open("search.log", "wt")


# Iterate over reads in BAM
for s in infile:
    t = robust_get_tag(s, "CB")          # cell barcode (if present)
    read_name = s.query_name
    align_chr = s.reference_name
    align_start = s.reference_start
    align_seq = s.query_sequence

    mystart = str(snv_info[lower_index]["pos"])
    read_tol += 1
    read_len = len(align_seq)

    # Progress logging
    if read_tol % 1000000 == 0:
        print("scanning read " + str(read_tol))

    lock = 0

    # Iterate through SNVs that may overlap this read
    for i in range(lower_index, snv_tol, 1):
        snv_pos = snv_info[i]["pos"]

        if snv_pos >= align_start:
            # Update lower bound index to avoid rescanning earlier SNVs
            if lock == 0:
                lower_index = i
            lock += 1

            # Check if SNV lies within read span
            if snv_pos <= align_start + read_len:
                read_cover += 1

                # Reference motif
                motif_pos = snv_info[i]["motif_pos"]

                # Alternative motif (ref base replaced by alt base)
                motif_neg = (
                    motif_pos[:motif_len]
                    + snv_info[i]["alt_allele"]
                    + motif_pos[motif_len + 1:]
                )

                # Direct full-length motif match
                if re.search(motif_pos, align_seq):
                    read_wild += 1
                elif re.search(motif_neg, align_seq):
                    read_mutated += 1
                else:
                    # Handle edge cases where motif overlaps read boundary
                    delta = snv_pos - align_start

                    # SNV near right end of read
                    if read_len - delta < motif_len:
                        motif_pos_part = motif_pos[0:motif_len + 1]
                        motif_neg_part = motif_neg[0:motif_len + 1]
                        seq_part = align_seq[read_len - 2 * motif_len - 1:read_len]

                        if re.search(motif_pos_part, seq_part):
                            read_wild += 1
                        elif re.search(motif_neg_part, seq_part):
                            read_mutated += 1

                    # SNV near left end of read
                    elif delta <= motif_len:
                        motif_pos_part = motif_pos[motif_len:len(motif_pos)]
                        motif_neg_part = motif_neg[motif_len:len(motif_neg)]
                        seq_part = align_seq[0:2 * motif_len + 1]

                        if re.search(motif_pos_part, seq_part):
                            read_wild += 1
                        elif re.search(motif_neg_part, seq_part):
                            read_mutated += 1

                    # Ambiguous case: no allele confidently detected
                    else:
                        fp.write(str(snv_info[i]) + "\n")
                        fp.write(
                            str(align_chr) + ":" +
                            str(align_start) + ":" +
                            align_seq + "\n"
                        )
                        read_noAllele += 1
            else:
                break

infile.close()
fp.close()

# Print summary statistics
print(
    str(read_tol) + ":" +
    str(read_cover) + ":" +
    str(read_wild) + ":" +
    str(read_mutated) + ":" +
    str(read_noAllele)
)
