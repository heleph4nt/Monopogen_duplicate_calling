#!/usr/bin/env python3
"""
The main interface of scPopGene
"""
# This file provides core BAM-level preprocessing utilities used in the
# Monopogen / scPopGene pipeline.
# Its responsibility is to normalize BAM headers, split BAMs by single cells,
# and perform pooled variant calling as an initial SNV discovery step.

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
from pysam import VariantFile

# Path to internal pipeline libraries
LIB_PATH = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "pipelines/lib"))

# Ensure internal libraries are discoverable
if LIB_PATH not in sys.path:
    sys.path.insert(0, LIB_PATH)

# Base directories for the pipeline
PIPELINE_BASEDIR = os.path.dirname(os.path.realpath(sys.argv[0]))
CFG_DIR = os.path.join(PIPELINE_BASEDIR, "cfg")

# Global logger for debugging and traceability
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)



def addChr(args):
    # This function standardizes chromosome names in a BAM file by
    # adding the "chr" prefix (e.g., 1 → chr1).
    # This prevents downstream incompatibilities between BAM files,
    # reference genomes, and tools such as samtools or bcftools.

    in_bam = args.bamFile
    prefix = 'chr'
    out_bam = in_bam + "tmp.bam"

    input_bam = pysam.AlignmentFile(in_bam, "rb")

    # Extract and modify BAM header
    new_head = input_bam.header.to_dict()
    for seq in new_head['SQ']:
        seq['SN'] = prefix + seq['SN']

    # Write a new BAM file with updated header and reads
    with pysam.AlignmentFile(out_bam, "wb", header=new_head) as outf:
        for read in input_bam.fetch():
            prefixed_chrom = prefix + read.reference_name

            # Reconstruct each aligned read explicitly
            a = pysam.AlignedSegment(outf.header)
            a.query_name = read.query_name
            a.query_sequence = read.query_sequence
            a.reference_name = prefixed_chrom
            a.flag = read.flag
            a.reference_start = read.reference_start
            a.mapping_quality = read.mapping_quality
            a.cigar = read.cigar
            a.next_reference_id = read.next_reference_id
            a.next_reference_start = read.next_reference_start
            a.template_length = read.template_length
            a.query_qualities = read.query_qualities
            a.tags = read.tags

            outf.write(a)

    input_bam.close()
    outf.close()

    # Index the new BAM and replace the original
    os.system("samtools index " + out_bam)
    os.system("mv " + out_bam + " " + in_bam)
    os.system("mv " + out_bam + ".bai " + in_bam + ".bai")




def sort_chr(chr_lst):
    # Sort chromosome identifiers in biological order (1–22),
    # supporting both numeric ("1") and prefixed ("chr1") formats.

    chr_lst_sort = []
    for i in range(1, 23):
        i = str(i)
        if i in chr_lst:
            chr_lst_sort.append(i)
        i_chr = "chr" + i
        if i_chr in chr_lst:
            chr_lst_sort.append(i_chr)

    return chr_lst_sort




def bamSplit(para):
    # Split a merged BAM file into per-cell BAM files.
    # This is a core step for single-cell processing, ensuring that
    # each BAM represents exactly one cell with a clean read group.

    para_lst = para.strip().split(":")
    chr = para_lst[0]      # chromosome (not directly used here)
    cell = para_lst[1]     # cell barcode / cell ID
    out = para_lst[2]      # output directory
    app_path = para_lst[3]

    samtools = app_path + "/samtools"
    output_bam = out + "/Bam/merge.filter.targeted.bam"

    infile = pysam.AlignmentFile(output_bam, "rb")

    # Modify read group (RG) to represent a single cell
    tp = infile.header.to_dict()
    if len(tp['RG']) > 1:
        tp['RG'] = [tp['RG'][0]]

    tp['RG'][0]['SM'] = cell
    tp['RG'][0]['ID'] = cell

    cnt = 0
    outfile = pysam.AlignmentFile(
        out + "/Bam/split_bam/" + cell + ".bam", "wb", header=tp)

    for s in infile:
        # Attempt to retrieve the cellular barcode from CB tag
        t = robust_get_tag(s, "CB")

        if t != "NotFound":
            # Standard case: explicit cell barcode
            if t == cell:
                outfile.write(s)
                cnt += 1
        else:
            # Fallback: cell ID embedded in read name
            if re.search(cell, s.query_name):
                outfile.write(s)
                cnt += 1

    outfile.close()
    infile.close()

    # Index the per-cell BAM
    os.system(samtools + " index " + out + "/Bam/split_bam/" + cell + ".bam")

    return cnt




def jointCall(para):
    # Perform joint SNV calling across multiple single-cell BAMs.
    # This corresponds to the pooled pileup SNV discovery step
    # described in the Monopogen workflow.

    para_lst = para.strip().split(">")
    jobid = para_lst[0]     # genomic region or chunk ID
    chr = para_lst[1]
    out = para_lst[2]
    app_path = para_lst[3]
    reference = para_lst[4]

    samtools = app_path + "/samtools"
    bcftools = app_path + "/bcftools"
    bgzip = app_path + "/bgzip"

    # List of per-cell BAM files
    bam_filter = out + "/Bam/split_bam/cell.bam.lst"

    # Construct mpileup + normalization command
    cmd1 = samtools + " mpileup -b " + bam_filter + \
           " -f " + reference + \
           " -r " + jobid + \
           " -q 20 -Q 20 -t DP4 -d 10000 -v "

    cmd1 += " | " + bcftools + " view "
    cmd1 += " | " + bcftools + " norm -m-both -f " + reference
    cmd1 += " | grep -v \"<X>\" | grep -v INDEL | "
    cmd1 += bgzip + " -c > " + out + "/somatic/" + jobid + ".cell.gl.vcf.gz"

    print(cmd1)
    output = os.system(cmd1)

    # Optional cleanup of temporary BAM files (currently disabled)
    f = open(bam_filter, "r")
    for x in f:
        x = x.strip()
        # os.system("rm " + x)
        # os.system("rm " + x + ".bai")
    f.close()

    if output == 0:
        return jobid
