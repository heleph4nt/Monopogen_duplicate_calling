#!/usr/bin/env python3

# This file implements the germline SNV calling pipeline in Monopogen.
# It handles input validation, dependency checking, BAM preprocessing,
# and filtering of reads prior to LD-based genotyping and refinement.

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
from bamProcess import *
import multiprocessing as mp
from multiprocessing import Pool


# Path to internal pipeline libraries
LIB_PATH = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "pipelines/lib"))

# Ensure internal libraries are discoverable
if LIB_PATH not in sys.path:
    sys.path.insert(0, LIB_PATH)

# Pipeline directories
PIPELINE_BASEDIR = os.path.dirname(os.path.realpath(sys.argv[0]))
CFG_DIR = os.path.join(PIPELINE_BASEDIR, "cfg")


# Global logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)


def print_parameters_given(args):
    # Print all user-specified parameters (except function handler)
    # Useful for reproducibility and debugging.
    logger.info("Parameters in effect:")
    for arg in vars(args):
        if arg == "func":
            continue
        logger.info("--{} = [{}]".format(arg, vars(args)[arg]))


def validate_sample_list_file(args):
    # Validate the sample list file, which specifies:
    # sample_name, absolute_path_to_bam, contamination_rate
    #
    # This function ensures:
    # - BAM files exist and are indexed
    # - Absolute paths are used
    # - Contamination rates are valid
    # - Optional check for improper hard-clipped reads

    if args.check_hard_clipped:
        out = os.popen("command -v bioawk").read().strip()
        assert out != "", "Program bioawk cannot be found!"

    assert os.path.isfile(args.sample_list), \
        "Sample index file {} cannot be found!".format(args.sample_list)

    try:
        with open(args.sample_list) as f_in:
            for line in f_in:
                record = line.strip().split("\t")
                logger.debug("Checking sample {}".format(record[0]))

                # Each line must have exactly 3 fields
                assert len(record) == 3, \
                    "Sample {} does not have exactly 3 columns!".format(record[0])

                # BAM file existence and indexing
                assert os.path.isfile(record[1]), \
                    "Bam file {} cannot be found!".format(record[1])
                assert os.path.isfile(record[1] + ".bai"), \
                    "Bam file {} has not been indexed!".format(record[1])
                assert os.path.isabs(record[1]), \
                    "Please use absolute path for bam file {}!".format(record[1])

                # Optional check for improper hard-clipped reads
                if args.check_hard_clipped:
                    logger.debug("Checking existence of hard-clipped reads.")
                    cmd = (
                        "samtools view {} | bioawk -c sam "
                        "'BEGIN {{count=0}} ($cigar ~ /H/)&&(!and($flag,256)) "
                        "{{count++}} END {{print count}}'"
                    ).format(record[1])
                    out = os.popen(cmd).read().strip()
                    assert out == "0", \
                        "Bam file {} contains hard-clipped reads!".format(record[1])

                # Validate contamination rate
                try:
                    cr = float(record[2])
                    assert 0.0 <= cr <= 1.0
                except:
                    logger.error(
                        "Invalid contamination rate for sample {}: {}".format(
                            record[0], record[2]
                        )
                    )
                    exit(1)

    except Exception:
        logger.error("There is something wrong with the sample index file.")
        print(sys.exc_info())
        raise sys.exc_info()[0]


def validate_user_setting_germline(args):
    # Validate user inputs required for germline SNV calling:
    # - Reference genome
    # - Imputation (LD) panel
    # - Region file
    # - Per-chromosome BAM lists

    assert os.path.isfile(args.reference), \
        "Reference genome {} cannot be found!".format(args.reference)
    assert os.path.isdir(args.imputation_panel), \
        "Imputation panel {} cannot be found!".format(args.imputation_panel)
    assert os.path.isfile(args.region), \
        "Region file {} cannot be found!".format(args.region)

    # Check per-chromosome BAM lists
    for chr in range(1, 23):
        bamFile = args.out + "/Bam/chr" + str(chr) + ".filter.bam.lst"
        with open(bamFile) as f_in:
            for line in f_in:
                line = line.strip()
                assert os.path.isfile(line), \
                    "Bam file {} cannot be found!".format(line)
                assert os.path.isfile(line + ".bai"), \
                    "Index file {}.bai cannot be found!".format(line)

    # Validate region file format
    with open(args.region) as f_in:
        for line in f_in:
            record = line.strip().split(",")
            assert len(record) == 3 or len(record) == 1, \
                "Invalid region specification: {}".format(line)


def check_dependencies(args):
    # Check availability of required external tools
    programs_to_check = (
        "vcftools", "bgzip", "bcftools",
        "beagle.27Jul16.86a.jar", "samtools",
        "picard.jar", "java"
    )

    for prog in programs_to_check:
        out = os.popen("command -v {}".format(args.app_path + "/" + prog)).read()
        assert out != "", "Program {} cannot be found!".format(prog)


def addChr(in_bam, samtools):
    # Add "chr" prefix to contig names in a BAM file.
    # This ensures consistency with reference genomes and panels.

    prefix = 'chr'
    out_bam = in_bam + "tmp.bam"

    input_bam = pysam.AlignmentFile(in_bam, "rb")
    new_head = input_bam.header.to_dict()

    for seq in new_head['SQ']:
        seq['SN'] = prefix + seq['SN']

    with pysam.AlignmentFile(out_bam, "wb", header=new_head) as outf:
        for read in input_bam.fetch():
            a = pysam.AlignedSegment(outf.header)
            a.query_name = read.query_name
            a.query_sequence = read.query_sequence
            a.reference_name = prefix + read.reference_name
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

    os.system(samtools + " index " + out_bam)
    os.system("mv " + out_bam + " " + in_bam)
    os.system("mv " + out_bam + ".bai " + in_bam + ".bai")


def BamFilter(myargs):
    # Filter reads from a BAM file for a given chromosome based on
    # maximum allowed mismatches (NM / nM tag).
    #
    # This reduces sequencing errors prior to SNV calling.

    bamFile = myargs.get("bamFile")
    search_chr = myargs.get("chr")
    samtools = myargs.get("samtools")
    id = myargs.get("id")
    max_mismatch = myargs.get("max_mismatch")
    out = myargs.get("out")

    os.system("mkdir -p " + out + "/Bam")

    infile = pysam.AlignmentFile(bamFile, "rb")

    # Detect whether contigs already have "chr" prefix
    contig_names = infile.references
    has_chr_prefix = any(contig.startswith("chr") for contig in contig_names)

    if not has_chr_prefix:
        logger.info(
            "Contig {} does not contain 'chr' prefix; adjusting.".format(search_chr)
        )
        search_chr = search_chr[3:]

    # Ensure proper read group (RG) information
    tp = infile.header.to_dict()
    sampleID = os.path.splitext(os.path.basename(bamFile))[0]
    tp.update({'RG': [{
        'SM': sampleID,
        'ID': sampleID,
        'LB': "0.1",
        'PL': "ILLUMINA",
        'PU': sampleID
    }]})

    outfile = pysam.AlignmentFile(
        out + "/Bam/" + id + "_" + search_chr + ".filter.bam",
        "wb",
        header=tp
    )

    for s in infile.fetch(search_chr):
        if s.has_tag("NM"): # Edit distance tag: The Levenstein distance between read and reference: which means how many bases are different between read and reference
            val = s.get_tag("NM")
        elif s.has_tag("nM"):
            val = s.get_tag("nM")
        else:
            continue # we ignore all of the reads that don't have any information about Levenstein distance.

        if val < max_mismatch: # if the distance is too large, we ignore it.
            outfile.write(s) 

    infile.close()
    outfile.close()

    os.system(samtools + " index " +
              out + "/Bam/" + id + "_" + search_chr + ".filter.bam")

    # Add "chr" prefix if needed
    if not has_chr_prefix:
        addChr(out + "/Bam/" + id + "_" + search_chr + ".filter.bam", samtools)

    return out + "/Bam/" + id + "_" + search_chr + ".filter.bam"


def robust_get_tag(read, tag_name):
    # Safe retrieval of BAM tags
    try:
        return read.get_tag(tag_name)
    except KeyError:
        return "NotFound"


def runCMD(cmd):
    # Execute a shell command
    output = os.system(cmd)
    if output == 0:
        return region