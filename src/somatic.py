#!/usr/bin/env python3

# Somatic SNV detection and LD refinement module of Monopogen

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


# Count number of SNVs in a VCF file (used for QC / diagnostics)
def withSNVs(invcf, path):
	print(invcf)
	os.system(path + "/tabix -p vcf " + invcf)
	vcf = VariantFile(invcf) 
	cnt = 0
	for rec in vcf.fetch():
		cnt += 1
	return cnt


# Execute a shell command
def runCMD(cmd):
	print(cmd)
	os.system(cmd)


# Safely retrieve a BAM tag (e.g. cell barcode)
def robust_get_tag(read, tag_name):  
	try:  
		return read.get_tag(tag_name)
	except KeyError:
		return "NotFound"


# Validate inputs for somatic variant calling
def validate_user_setting_somatic(args):

	assert os.path.isdir(args.out), "Germline output folder not found"
	assert os.path.isfile(args.region), "Region file not found"
	
	with open(args.region) as f_in:
		for line in f_in:
			record = line.strip().split(",")
			if len(record) == 1:
				jobid = record[0]
			else:
				logger.error("Only whole chromosome somatic calling is allowed")
				exit(1)

			# Required germline outputs
			bam_filter = args.out + "/Bam/" + record[0] + ".filter.bam.lst"
			gl_vcf = args.out + "/germline/" + jobid + ".gl.vcf.gz"
			phased_vcf = args.out + "/germline/" + jobid + ".phased.vcf.gz"

			assert os.path.isfile(bam_filter)
			assert os.path.isfile(gl_vcf)
			assert os.path.isfile(phased_vcf)
			assert os.path.isfile(args.barcode)


# Extract INFO field robustly from a VCF record
def getInfo_robust(rec, info):
	info_dt = rec.info.get(info)
	if info_dt is not None:
		if not isinstance(info_dt, float):
			info_dt = info_dt[0]
		info_dt = round(info_dt, 2)
	return info_dt


# Extract reads overlapping high-confidence SNV regions
def bamExtract(para):

	chr, out, app_path = para.strip().split(">")
	samtools = os.path.abspath(app_path) + "/samtools" 
	out = os.path.abspath(out)

	inbam = getBamName(chr, out)
	outbam = out + "/Bam/" + chr + ".filter.targeted.bam"
	chr_bed = out + "/somatic/" + chr + ".gl.vcf.filter.hc.bed"

	cmd = samtools + " view -b -L " + chr_bed + " " + inbam + " -o " + outbam

	with open(out + "/Script/bamExtract_" + chr + ".sh", "w") as f:
		f.write(cmd + "\n")
		f.write(samtools + " index " + outbam + "\n")

	os.system("bash " + out + "/Script/bamExtract_" + chr + ".sh")


# Extract germline feature information for somatic filtering
def featureInfo(para):

	region, out, app = para.strip().split(">")

	# Index germline VCFs
	pysam.tabix_index(out + "/germline/" + region + ".phased.vcf.gz", preset="vcf", force=True)
	pysam.tabix_index(out + "/germline/" + region + ".gl.vcf.gz", preset="vcf", force=True)

	# Load phased genotypes
	vcf_in = VariantFile(out + "/germline/" + region + ".phased.vcf.gz") 
	info_GT = {}

	for rec in vcf_in.fetch():
		GT = [value['GT'] for value in rec.samples.values()][0]
		id = f"{rec.chrom}:{rec.pos}:{rec.ref}:{rec.alts[0]}"
		info_GT[id] = GT

	# Prepare output files
	gl_vcf_in = VariantFile(out + "/germline/" + region + ".gl.vcf.gz") 
	gl_vcf_dp4 = open(out + "/somatic/" + region + ".gl.vcf.DP4", "w")
	gl_vcf_filter_dp4 = open(out + "/somatic/" + region + ".gl.vcf.filter.DP4", "w")
	gl_vcf_filter_bed = open(out + "/somatic/" + region + ".gl.vcf.filter.hc.bed", "w")
	gl_vcf_filter_txt = open(out + "/somatic/" + region + ".gl.vcf.filter.hc.pos", "w")

	germline_cutoff = 4
	depth_filter_novelSNV = 10

	# Loop over germline variants
	for rec in gl_vcf_in.fetch():
		info_I16 = rec.info.get('I16')
		info_QS = getInfo_robust(rec, 'QS')
		info_VDB = getInfo_robust(rec, 'VDB')
		info_RPB = getInfo_robust(rec, 'RPB')
		info_MQB = getInfo_robust(rec, 'MQB')  
		info_BQB = getInfo_robust(rec, 'BQB') 
		info_MQSB = getInfo_robust(rec, 'MQSB')
		info_SGB = getInfo_robust(rec, 'SGB')
		info_MQ0F = getInfo_robust(rec, 'MQ0F')

		id = f"{rec.chrom}:{rec.pos}:{rec.ref}:{rec.alts[0]}"
		gt_info_var = "NA"
		if id in info_GT:
			gt_info_var = f"{info_GT[id][0]}|{info_GT[id][1]}"

		# Write DP4-style summary
		line = f"{rec.chrom}\t{rec.pos}\t{rec.ref}\t{rec.alts[0]}\t{rec.info['DP']}\t" \
		       f"{info_I16[0]}\t{info_I16[1]}\t{info_I16[2]}\t{info_I16[3]}\t{gt_info_var}\t" \
		       f"{info_VDB}\t{info_QS}\t{info_RPB}\t{info_MQB}\t{info_BQB}\t{info_MQSB}\t{info_SGB}\t{info_MQ0F}\n"
		gl_vcf_dp4.write(line)

		ref_depth = info_I16[0] + info_I16[1]
		alt_depth = info_I16[2] + info_I16[3]

		# High-confidence filtering for somatic analysis
		if ((ref_depth >= 4 and alt_depth >= 4) or (id in info_GT)):
			gl_vcf_filter_bed.write(f"{rec.chrom}\t{rec.pos-1}\t{rec.pos}\n")
			gl_vcf_filter_txt.write(f"{rec.chrom}\t{rec.pos}\n")
			gl_vcf_filter_dp4.write(line)

	gl_vcf_dp4.close()
	gl_vcf_filter_dp4.close()
	gl_vcf_filter_txt.close()
	gl_vcf_filter_bed.close()

	# Extract reads overlapping somatic candidate sites
	bamExtract(para)
	return region
