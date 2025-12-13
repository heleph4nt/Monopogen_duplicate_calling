#!/usr/bin/env python
"""
This script implements the MonoVar probabilistic model for SNV calling
from single-cell sequencing data.

It reads pileup-like input from stdin and outputs a VCF file with
genotype likelihoods and variant calls across multiple single cells.
"""

import sys
import re
import pysam
import time
import fileinput
import copy
import json
import _pickle as cPickle
import math
import copyreg as copy_reg
import types
import multiprocessing as mp
from functools import partial
from operator import add
from contextlib import closing

# Import core MonoVar / Monopogen components
from utils import Utils_Functions
from alleles_prior import allele_prior
from Single_Cell_Ftrs_Pos import Single_Cell_Ftrs_Pos
from nu_prob_mat import Prob_matrix
from calc_variant_prob import Calc_Var_Prob
from genotype_prob_mat import Genotype_Prob_matrix
from nu_genotype_single_cell import Single_cell_genotype_records
from mp_genotype import MP_single_cell_genotype
from hzvcf import VCFRecord, VCFDocument, VRecord, VCF


# Required to allow multiprocessing to pickle class methods
def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)


# Instantiate utility and multiprocessing genotype helper
U = Utils_Functions()
M = MP_single_cell_genotype()


# -----------------------------
# Default parameter values
# -----------------------------
pe = 0.002          # sequencing error probability
pad = 0.2           # allelic dropout probability
thr = 0.05          # SNV calling threshold
m_thread = 1        # number of threads
CF_flag = 1         # consensus filter flag

input_args = {}


# -----------------------------
# Parse command-line arguments
# -----------------------------
argc = len(sys.argv)
i = 1
while (i < argc):
    if (sys.argv[i] == '-n'):
        n_cells = int(sys.argv[i + 1])
    elif (sys.argv[i] == '-p'):
        pe = float(sys.argv[i + 1])
    elif (sys.argv[i] == '-d'):
        pd = float(sys.argv[i + 1])     # deamination error (not always used)
    elif (sys.argv[i] == '-a'):
        pad = float(sys.argv[i + 1])
    elif (sys.argv[i] == '-f'):
        ref_file = sys.argv[i + 1]
        input_args['-f'] = 'Provided'
    elif (sys.argv[i] == '-b'):
        bam_file_list = sys.argv[i + 1]
        input_args['-b'] = 'Provided'
    elif (sys.argv[i] == '-o'):
        outfile = sys.argv[i + 1]
        input_args['-o'] = 'Provided'
    elif (sys.argv[i] == '-t'):
        thr = float(sys.argv[i + 1])
    elif (sys.argv[i] == '-c'):
        CF_flag = int(sys.argv[i + 1])
    elif (sys.argv[i] == '-m'):
        m_thread = int(sys.argv[i + 1])
    i = i + 2


# -----------------------------
# Mandatory argument checks
# -----------------------------
if '-f' not in input_args:
    print("Error: Reference genome file not provided.")
    exit(3)

if '-b' not in input_args:
    print("Error: List of BAM files not provided.")
    exit(3)

if '-o' not in input_args:
    print("Error: Output file not provided.")
    exit(3)

if CF_flag > 1:
    print("CF_flag must be 0 or 1.")
    exit(3)


# -----------------------------
# Extract cell/sample IDs from BAM files
# -----------------------------
bam_id_list = []
with open(bam_file_list) as f_bam_list:
    for filename in f_bam_list:
        filename = filename.strip()
        bam_id = U.Get_BAM_RG(filename)
        bam_id_list.append(bam_id)

n_cells = len(bam_id_list)


# -----------------------------
# Initialize multiprocessing
# -----------------------------
pool = mp.Pool(processes=m_thread)


# -----------------------------
# Global constants
# -----------------------------
cell_no_threshold = n_cells / 2
max_allele_cnt = 2 * n_cells + 1     # possible total alt allele counts
theta = 0.001                        # heterozygosity rate
Base_dict = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}


# -----------------------------
# Precompute combinatorics
# -----------------------------
factorial_list = U.Create_Factorial_List(max_allele_cnt)
nCr_matrix = U.Create_nCr_mat(max_allele_cnt, factorial_list)


# -----------------------------
# Prior on number of mutated cells
# -----------------------------
prior_variant_dict = {}
for i in range(n_cells + 1):
    prior_variant_dict[i] = U.calc_prior(theta, i, 1)


# -----------------------------
# Initialize VCF output
# -----------------------------
f_vcf = open(outfile, 'w')
vcf = VCFDocument(f_vcf)
vcf.populate_fields(bam_id_list)
vcf.populate_reference(ref_file)
vcf.print_header()


# -----------------------------
# Allocate reusable structures
# -----------------------------
All_single_cell_ftrs_list = n_cells * [None]
read_flag_row = n_cells * [None]
alt_allele_flag_row = n_cells * [None]


# ==========================================================
# Main loop: one genomic position per input line
# ==========================================================
for line in sys.stdin:
    row = line.strip().split('\t')
    contig = row[0]
    pos = int(row[1])
    refBase = U.refineBase(row[2])

    total_depth = 0
    total_ref_depth = 0

    # Build Single_Cell_Ftrs_Pos object for each cell
    for i in range(1, n_cells + 1):
        curr_cell = Single_Cell_Ftrs_Pos(refBase, row[3 * i: 3 * i + 3])
        total_depth += curr_cell.depth
        total_ref_depth += curr_cell.refDepth
        All_single_cell_ftrs_list[i - 1] = curr_cell

    Alt_count = total_depth - total_ref_depth
    Alt_freq = 0 if total_depth == 0 else float(Alt_count) / total_depth


    # -----------------------------
    # Pre-filters to remove noise
    # -----------------------------
    if total_ref_depth == total_depth:
        continue
    elif (total_depth > 30) and ((Alt_count <= 2) or (Alt_freq <= 0.001)):
        continue
    elif Alt_freq <= 0.01:
        continue
    elif refBase not in ['A', 'T', 'G', 'C']:
        continue
    elif total_depth <= 10:
        continue


    # -----------------------------
    # Collect read-supported cells
    # -----------------------------
    read_supported_cell_list = []
    total_alt_allele_count = [0, 0, 0, 0]
    cell_index_list = []
    info_list = ['GT:AD:DP:GQ:PL']

    c = 1
    for j in range(n_cells):
        cell = All_single_cell_ftrs_list[j]
        read_flag = U.checkReadPresence(cell)

        if read_flag:
            alt_flag = U.CheckAltAllele(cell)
            cell.Get_base_call_string_nd_quals(refBase)
            cell.Get_Alt_Allele_Count()
            cell.Set_Cell_Index(c)

            if alt_flag:
                total_alt_allele_count = U.update_alt_count(
                    total_alt_allele_count, cell)

            read_supported_cell_list.append(cell)
            cell_index_list.append(c)
        else:
            alt_flag = 0

        read_flag_row[j] = read_flag
        alt_allele_flag_row[j] = alt_flag
        c += 1


    read_smpl_count = sum(read_flag_row)
    alt_smpl_count = sum(alt_allele_flag_row)
    Alt_count = max(total_alt_allele_count)

    if Alt_count == 0:
        continue

    altBase = Base_dict[total_alt_allele_count.index(Alt_count)]


    # -----------------------------
    # Prior allele matrix
    # -----------------------------
    prior_allele_mat = U.Get_prior_allele_mat(
        read_smpl_count, alt_smpl_count,
        cell_no_threshold, total_depth, Alt_freq, pe)


    for cell in read_supported_cell_list:
        cell.Store_Addl_Info(refBase, altBase, Alt_freq, prior_allele_mat)


    # -----------------------------
    # Compute P(no variant)
    # -----------------------------
    prior_variant_number = prior_variant_dict[len(read_supported_cell_list)]
    Calc_var_prob_obj = Calc_Var_Prob(
        read_supported_cell_list, prior_allele_mat, len(read_supported_cell_list))

    (zero_variant_prob, denominator) = Calc_var_prob_obj.Calc_Zero_Var_Prob(
        n_cells, 10000, nCr_matrix, pad, prior_variant_number)


    # -----------------------------
    # Variant decision
    # -----------------------------
    if zero_variant_prob <= thr:

        (max_prob_ratio, _) = U.find_max_prob_ratio(
            Calc_var_prob_obj.matrix, Calc_var_prob_obj.matrix_shape)

        (oddsRatio, _) = U.calc_strand_bias(
            read_supported_cell_list, Alt_count)

        baseQranksum = U.Calc_Base_Q_Rank_Sum(
            read_supported_cell_list, refBase, altBase)[0]

        genotype_dict = {0: '0/0', 1: '0/1', 2: '1/1'}

        # Parallel per-cell genotype inference
        func = partial(
            M.get_info_string,
            read_supported_cell_list,
            prior_allele_mat,
            n_cells,
            nCr_matrix,
            prior_variant_number,
            denominator,
            genotype_dict)

        output = pool.map(func, range(len(read_supported_cell_list)))

        # Build VCF barcode and FORMAT fields
        barcode = '<'
        read_supported_info_list = [p[0] for p in output]
        read_supported_barcodes = [p[1] for p in output]

        for j in range(n_cells):
            if All_single_cell_ftrs_list[j].depth == 0:
                info_list.append('./.')
                barcode += 'X'
            else:
                info_list.append(read_supported_info_list.pop(0))
                barcode += read_supported_barcodes.pop(0)

        barcode += '>'

        Qual = 323 if zero_variant_prob == 0 else -math.log10(zero_variant_prob)

        (AC, AF, AN) = U.Calc_chr_count(barcode)
        QD = U.Calc_qual_depth(barcode, All_single_cell_ftrs_list, Qual)
        PSARR = U.Calc_Per_Smpl_Alt_Ref_Ratio(
            total_ref_depth, Alt_count, read_smpl_count, alt_smpl_count)

        vcf_record = VRecord(contig, pos)
        info_names = [i[0] for i in vcf.info_fields]
        info_values = [
            AC, "%.2f" % AF, AN,
            "%.2f" % baseQranksum,
            total_depth, "%.2f" % QD,
            "%.2f" % oddsRatio,
            "%.2f" % max_prob_ratio,
            "%.2f" % PSARR
        ]

        info_string = ';'.join('%s=%s' % t for t in zip(info_names, info_values))

        # Consensus filtering
        if CF_flag == 1 and U.Consensus_Filter(barcode):
            vcf_record.get6fields(refBase, altBase, '.', Qual, 'PASS', info_string)
        else:
            vcf_record.get6fields(refBase, altBase, '.', Qual, '.', info_string)

        vcf_record.format_vcf(info_list)
        vcf_record.get_passcode(barcode)
        vcf.print_my_record(vcf_record)
