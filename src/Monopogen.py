#!/usr/bin/env python3
"""
The main interface of monopogen
"""

# This file is the main entry point of Monopogen.
# It defines the command-line interface and orchestrates the full pipeline:
# preprocessing → germline calling → somatic calling.

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

# Import pipeline modules
from bamProcess import *
from germline import *
from somatic import *

import multiprocessing as mp
from multiprocessing import Pool


# Internal library path
LIB_PATH = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "pipelines/lib"))

if LIB_PATH not in sys.path:
    sys.path.insert(0, LIB_PATH)

PIPELINE_BASEDIR = os.path.dirname(os.path.realpath(sys.argv[0]))
CFG_DIR = os.path.join(PIPELINE_BASEDIR, "cfg")


# Global logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    '[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)



def error_check(all, output, step):
    # Check whether all expected jobs completed successfully.
    # If any job ID is missing from the output list, abort execution.
    job_fail = 0
    for id in all:
        if id not in output:
            logger.error("In " + step + " step " + id + " failed!")
            job_fail += 1

    if job_fail > 0:
        logger.error("Failed! See instructions above.")
        exit(1)



def germline(args):
    # Run the germline variant calling pipeline.
    # This includes:
    # - mpileup-based SNV scanning
    # - genotype likelihood computation
    # - LD-based refinement using a population panel (e.g. 1KG)

    logger.info("Performing germline variant calling...")
    print_parameters_given(args)

    logger.info("Checking existence of essential resource files...")
    validate_user_setting_germline(args)

    logger.info("Checking dependencies...")
    check_dependencies(args)

    out = args.out
    os.system("mkdir -p " + out)
    os.system("mkdir -p " + out + "/germline")
    os.system("mkdir -p " + out + "/Script")

    joblst = []

    # Iterate over genomic regions to be processed
    with open(args.region) as f_in:
        for line in f_in:
            record = line.strip().split(",")

            # Define region/job ID
            if len(record) == 1:
                jobid = record[0]
            elif len(record) == 3:
                jobid = record[0] + ":" + record[1] + "-" + record[2]

            # BAM list for this chromosome
            bam_filter = args.out + "/Bam/" + record[0] + ".filter.bam.lst"

            # Count number of samples
            N_sample = 0
            with open(bam_filter) as p:
                for s in p:
                    N_sample += 1

            # Population reference panel (1KG)
            imputation_vcf = (
                args.imputation_panel +
                "CCDG_14151_B01_GRM_WGS_2020-08-05_" +
                record[0] +
                ".filtered.shapeit2-duohmm-phased.vcf.gz"
            )

            # Step 1: pooled mpileup to generate genotype likelihoods
            cmd1 = (
                samtools + " mpileup -b " + bam_filter +
                " -f " + args.reference +
                " -r " + jobid +
                " -q 20 -Q 20 --incl-flags 0 --excl-flags 0 " +
                "-t DP -d 10000000 -v "
            )
            cmd1 += (
                " | " + bcftools + " view | " +
                bcftools + " filter -e 'REF !~ \"^[ATGC]$\"' | " +
                bcftools + " norm -m-both -f " + args.reference
            )
            cmd1 += (
                " | grep -v \"<X>\" | " + bgzip +
                " -c > " + args.out + "/germline/" +
                jobid + ".gl.vcf.gz"
            )

            # Step 2: LD-based genotype refinement using Beagle
            cmd3 = (
                java + " -Xmx20g -jar " + beagle +
                " gl=" + out + "/germline/" + jobid + ".gl.vcf.gz" +
                " ref=" + imputation_vcf +
                " chrom=" + record[0] +
                " out=" + out + "/germline/" + jobid + ".gp " +
                "impute=false modelscale=2 gprobs=true niterations=0"
            )

            # Step 3: phasing
            cmd5 = (
                java + " -Xmx20g -jar " + beagle +
                " gt=" + out + "/germline/" + jobid + ".germline.vcf" +
                " ref=" + imputation_vcf +
                " chrom=" + record[0] +
                " out=" + out + "/germline/" + jobid + ".phased " +
                "impute=false modelscale=2 gprobs=true niterations=0"
            )
            cmd5 += "\nrm " + out + "/germline/" + jobid + ".germline.vcf"

            # Write job script
            f_out = open(out + "/Script/runGermline_" + jobid + ".sh", "w")

            if args.step in ("varScan", "all"):
                f_out.write(cmd1 + "\n")

            if args.step in ("varImpute", "all"):
                f_out.write(cmd3 + "\n")
                cmd4 = (
                    "zless -S " + out + "/germline/" +
                    jobid + ".gp.vcf.gz > " +
                    out + "/germline/" + jobid + ".germline.vcf"
                )
                f_out.write(cmd4 + "\n")

            if args.step in ("varPhasing", "all"):
                f_out.write(cmd5 + "\n")

            f_out.close()
            joblst.append("bash " + out + "/Script/runGermline_" + jobid + ".sh")

    # Run jobs in parallel unless --norun is set
    if args.norun != "TRUE":
        with Pool(processes=args.nthreads) as pool:
            result = pool.map(runCMD, joblst)



def somatic(args):
    # Run somatic SNV calling using Monopogen’s LD-based refinement
    # and single-cell evidence aggregation.

    validate_user_setting_somatic(args)
    os.system("mkdir -p " + args.out + "/somatic")

    chr_lst = []
    region_lst = []

    with open(args.region) as f_in:
        for line in f_in:
            record = line.strip().split(",")
            if len(record) == 1:
                region = record[0]
            elif len(record) == 3:
                region = record[0] + ":" + record[1] + "-" + record[2]
            chr_lst.append(record[0])
            region_lst.append(region)

    chr_lst = list(set(chr_lst))

    # Step 1: extract feature-level information
    if args.step in ("featureInfo", "all"):
        logger.info("Get feature information from sequencing data...")
        joblst = [id + ">" + args.out + ">" + args.app_path for id in region_lst]
        with Pool(processes=args.nthreads) as pool:
            result = pool.map(featureInfo, joblst)
        error_check(all=chr_lst, output=result, step="featureInfo")

    # Step 2: extract single-cell read-level evidence
    if args.step in ("cellScan", "all"):
        logger.info("Collect single cell level information...")
        chr_lst = sort_chr(chr_lst)
        joblst = [
            id + ">" + args.out + ">" + args.reference + ">" + args.barcode
            for id in chr_lst
        ]
        with Pool(processes=args.nthreads) as pool:
            result = pool.map(bam2mat, joblst)
        error_check(all=chr_lst, output=result, step="cellScan")

    # Step 3: LD-based refinement
    if args.step in ("LDrefinement", "all"):
        logger.info("Run LD refinement...")
        joblst = [id + ">" + args.out + ">" + args.app_path for id in chr_lst]
        with Pool(processes=args.nthreads) as pool:
            result = pool.map(LDrefinement, joblst)
        error_check(all=chr_lst, output=result, step="LDrefinement")



def preProcess(args):
    # Preprocess BAM files by filtering reads with excessive mismatches
    # and splitting data by chromosome.

    logger.info("Performing data preprocess...")
    print_parameters_given(args)

    assert os.path.isfile(args.bamFile), \
        "The bam file {} cannot be found!".format(args.bamFile)

    out = args.out
    os.system("mkdir -p " + out)
    os.system("mkdir -p " + out + "/Bam")

    sample = []

    with open(args.bamFile) as f_in:
        for line in f_in:
            record = line.strip().split(",")
            sample.append(record[0])
            assert len(record) == 2
            assert os.path.isfile(record[1])
            assert os.path.isfile(record[1] + ".bai")

    para_lst = []

    with open(args.bamFile) as f_in:
        for line in f_in:
            record = line.strip().split(",")
            for chr in range(1, 23):
                para_single = dict(
                    chr="chr" + str(chr),
                    out=args.out,
                    id=record[0],
                    bamFile=record[1],
                    max_mismatch=args.max_mismatch,
                    samtools=samtools
                )
                para_lst.append(para_single)

    with Pool(processes=args.nthreads) as pool:
        pool.map(BamFilter, para_lst)

    # Write per-chromosome BAM lists
    for chr in range(1, 23):
        with open(args.out + "/Bam/chr" + str(chr) + ".filter.bam.lst", "w") as bamlist:
            for s in sample:
                bamlist.write(
                    args.out + "/Bam/" + s + "_chr" + str(chr) + ".filter.bam\n"
                )



def main():
    # Command-line interface definition

    parser = argparse.ArgumentParser(
        description="Monopogen: SNV calling from single cell sequencing",
        formatter_class=argparse.RawTextHelpFormatter
    )

    subparsers = parser.add_subparsers(dest="subcommand")

    # preProcess command
    parser_preProcess = subparsers.add_parser("preProcess")
    parser_preProcess.add_argument("-b", "--bamFile", required=True)
    parser_preProcess.add_argument("-o", "--out")
    parser_preProcess.add_argument("-a", "--app-path", required=True)
    parser_preProcess.add_argument("-m", "--max-mismatch", type=int, default=3)
    parser_preProcess.add_argument("-t", "--nthreads", type=int, default=1)
    parser_preProcess.set_defaults(func=preProcess)

    # germline command
    parser_germline = subparsers.add_parser("germline")
    parser_germline.add_argument("-r", "--region", required=True)
    parser_germline.add_argument(
        "-s", "--step", choices=["varScan", "varImpute", "varPhasing", "all"]
    )
    parser_germline.add_argument("-o", "--out")
    parser_germline.add_argument("-g", "--reference", required=True)
    parser_germline.add_argument("-p", "--imputation-panel", required=True)
    parser_germline.add_argument("-a", "--app-path", required=True)
    parser_germline.add_argument("-t", "--nthreads", type=int, default=1)
    parser_germline.add_argument("-n", "--norun", default="FALSE")
    parser_germline.set_defaults(func=germline)

    # somatic command
    parser_somatic = subparsers.add_parser("somatic")
    parser_somatic.add_argument("-i", "--input-folder", required=True)
    parser_somatic.add_argument("-r", "--region", required=True)
    parser_somatic.add_argument("-l", "--barcode", required=True)
    parser_somatic.add_argument("-a", "--app-path", required=True)
    parser_somatic.add_argument("-t", "--nthreads", type=int, default=22)
    parser_somatic.add_argument(
        "-s", "--step",
        choices=["featureInfo", "cellScan", "LDrefinement", "monovar", "all"]
    )
    parser_somatic.add_argument("-g", "--reference", required=True)
    parser_somatic.set_defaults(func=somatic)

    args = parser.parse_args()

    if args.subcommand is None:
        parser.print_help()
        exit(1)

    if args.subcommand == "somatic":
        args.out = args.input_folder

    # Global tool paths
    global out, samtools, bcftools, bgzip, java, beagle
    out = os.path.abspath(args.out)
    samtools = os.path.abspath(args.app_path) + "/samtools"
    bcftools = os.path.abspath(args.app_path) + "/bcftools"
    bgzip = os.path.abspath(args.app_path) + "/bgzip"
    java = "java"
    beagle = os.path.abspath(args.app_path) + "/beagle.27Jul16.86a.jar"

    args.func(args)

    logger.info("Success! See instructions above.")



if __name__ == "__main__":
    main()
