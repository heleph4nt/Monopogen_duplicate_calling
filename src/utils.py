"""
The MIT License
...
"""

import math
import heapq
import copy
import re
import sys
import numpy as np
import pysam
from operator import add
from scipy import stats
import fileinput
from contextlib import closing
from base_q_ascii import Ascii_Table
from alleles_prior import allele_prior

# ASCII table used to convert base quality characters into error probabilities
base_q_tbl = Ascii_Table()


class Utils_Functions:

    def nCr(self, n, r):
        # Compute binomial coefficient n choose r using factorials
        f = math.factorial
        return f(n) / f(r) / f(n - r)

    def find_second_smallest(self, l):
        # Return the second smallest element in a list
        return heapq.nsmallest(2, l)[-1]

    def GetPosition(self, row):
        # Extract genomic position from a parsed pileup row
        position = int(row[1])
        return position

    def GetRefBase(self, row):
        # Extract reference base from a parsed pileup row
        refbase = row[2]
        return refbase

    def single_cell_ftrs(self, row, cell_index):
        # Extract features for a specific cell from a pileup row
        newrow = row[cell_index + 2]
        newrow = newrow.split('\t')

        # First 6 entries are numeric features (depth, counts, etc.)
        cell_ftr = [int(i) for i in newrow[0:6]]

        # Read bases string (remove spaces)
        read_bases = newrow[6]
        read_bases = read_bases.replace(' ', '')

        # Base quality string
        newrow[7] = newrow[7].replace(' ', '')
        if newrow[7] == '[]':
            base_qual_list_f = []
        else:
            base_qual_list = newrow[7].split(',')
            base_qual_list[0] = base_qual_list[0].replace('[', '')
            base_qual_list[-1] = base_qual_list[-1].replace(']', '')
            base_qual_list_f = [float(i) for i in base_qual_list]

        # Return extracted features, base calls, and base qualities
        return [cell_ftr, read_bases, base_qual_list_f]

    def read_last_line(self, s):
        # Read last line of a file efficiently (used to get last genomic position)
        with open(s, "rb") as f:
            f.seek(-2, 2)
            while (f.read(1) != '\n'):
                f.seek(-2, 1)
            last = f.readline()
        last = last.replace('\n', '')
        pos = int(last.split('\t')[1])
        return pos

    def checkReadPresence(self, single_cell_dict):
        # Check whether a cell has any reads covering the position
        if single_cell_dict.depth == 0:
            return 0
        else:
            return 1

    def Create_Factorial_List(self, max_allele_cnt):
        # Precompute factorials up to max_allele_cnt
        factorial_list = max_allele_cnt * [None]
        f = math.factorial
        for i in range(max_allele_cnt):
            factorial_list[i] = f(i)
        return factorial_list

    def Create_nCr_mat(self, max_allele_cnt, factorial_list):
        # Precompute nCr values into a matrix for fast lookup
        ncr_mat = np.zeros((max_allele_cnt, max_allele_cnt))
        for i in range(max_allele_cnt):
            for j in range(max_allele_cnt):
                ncr_mat[j, i] = factorial_list[j] / \
                    (factorial_list[i] * factorial_list[j - i])
        return ncr_mat

    def CheckAltAllele(self, single_cell_dict):
        # Determine if the cell has any alternate allele observations
        alt_c = single_cell_dict.depth - single_cell_dict.refDepth
        if alt_c == 0:
            return 0
        else:
            return 1

    def RefCountString(self, read_base):
        # Count reference allele reads on forward (.) and reverse (,) strands
        forward_ref_c = read_base.count('.')
        reverse_ref_c = read_base.count(',')
        RefCount = forward_ref_c + reverse_ref_c
        return (forward_ref_c, reverse_ref_c, RefCount)

    def refineBase(self, b):
        # Normalize base character (remove spaces, uppercase)
        b = b.replace(' ', '')
        nb = b.upper()
        return nb

    def copy_list_but_one(self, i_list, index):
        # Return a copy of the list excluding one element (used for leave-one-out)
        nu_list = []
        len_i_list = len(i_list)
        len_nu_list = len_i_list - 1
        for i in range(index):
            nu_list.append(i_list[i])
        for i in range(index + 1, len_i_list):
            nu_list.append(i_list[i])
        return (nu_list, len_nu_list)

    def update_alt_count(self, alt_allele_count, cell_ftr_dict):
        # Update global alternate allele counts across cells
        alt_allele_count[0] += cell_ftr_dict.A_cnt
        alt_allele_count[1] += cell_ftr_dict.T_cnt
        alt_allele_count[2] += cell_ftr_dict.G_cnt
        alt_allele_count[3] += cell_ftr_dict.C_cnt
        return alt_allele_count

    def find_indel(self, string, pattern):
        # Find insertions or deletions using regex, remove them, and return cleaned string
        l = [x.group() for x in re.finditer(pattern, string)]
        len_indel = [int(re.split(r'(\d+)', i)[1]) for i in l]
        spans = [i.span() for i in re.finditer(pattern, string)]
        newspan = []
        for i in range(len(len_indel)):
            new_end = spans[i][0] + 1 + len_indel[i] + len(str(len_indel[i]))
            newspan.append((spans[i][0], new_end))
        final_indel_list = [string[i1:i2] for (i1, i2) in newspan]
        new_string = string
        for i in final_indel_list:
            new_string = new_string.replace(i, '')
        return (final_indel_list, new_string)

    def Alt_count(self, string):
        # Remove indels and count nucleotide occurrences
        cp_string = copy.deepcopy(string)
        (ins_list, ins_rmvd_str) = self.find_indel(
            cp_string, '\+[0-9]+[ACGTNacgtn]+')
        ins_count = len(ins_list)
        (del_list, del_ins_rmvd_str) = self.find_indel(
            ins_rmvd_str, '-[0-9]+[ACGTNacgtn]+')
        del_count = len(del_list)
        A_cnt = del_ins_rmvd_str.count('A') + del_ins_rmvd_str.count('a')
        T_cnt = del_ins_rmvd_str.count('T') + del_ins_rmvd_str.count('t')
        G_cnt = del_ins_rmvd_str.count('G') + del_ins_rmvd_str.count('g')
        C_cnt = del_ins_rmvd_str.count('C') + del_ins_rmvd_str.count('c')
        N_cnt = del_ins_rmvd_str.count('N') + del_ins_rmvd_str.count('n')
        return (del_ins_rmvd_str, ins_count, del_count, A_cnt, T_cnt, G_cnt, C_cnt, N_cnt)

    def Count_Start_and_End(self, s):
        # Count read starts (^) and ends ($) in pileup string
        end_counts = s.count('$')
        ns = s.replace('$', '')
        start_counts = 0
        i = 0
        fs = ''
        while (i < len(ns)):
            if ns[i] == '^':
                i += 2
                start_counts += 1
            else:
                fs = fs + ns[i]
                i += 1
        return (start_counts, end_counts, fs)

    def Create_base_call_string(self, s, ref):
        # Convert pileup symbols into explicit base calls
        l = ['.', ',', 'a', 'A', 'c', 'C', 't', 'T', 'g', 'G', '*']
        sn = ''
        for i in s:
            if i in l:
                sn = sn + i
        snn = ''
        for i in sn:
            if i == '.' or i == ',' or i == '*':
                snn = snn + ref
            elif i == 'a':
                snn = snn + 'A'
            elif i == 'c':
                snn = snn + 'C'
            elif i == 't':
                snn = snn + 'T'
            elif i == 'g':
                snn = snn + 'G'
            else:
                snn = snn + i
        return snn

    def Get_base_qual_list(self, s):
        # Convert base quality ASCII string into probability and integer lists
        len_s = len(s)
        base_q_list = len_s * [None]
        base_q_int_list = len_s * [None]
        for i in range(len_s):
            err_p = base_q_tbl.base_q_dict[s[i]]
            err_int = base_q_tbl.base_q_int_dict[s[i]]
            base_q_list[i] = err_p
            base_q_int_list[i] = err_int
        return (base_q_list, base_q_int_list)

    def find_min_list(self, actual_list, flag_list, last_min_loc, last_min_index):
        # Find minimum value with continuity constraint
        if ((sum(flag_list) == 1) & (actual_list[last_min_index] == last_min_loc + 1)):
            min_index = last_min_index
            min_loc = actual_list[last_min_index]
        else:
            min_loc = min(actual_list)
            min_index = actual_list.index(min_loc)
        return(min_loc, min_index)

    def feature_row(self, rlist):
        # Build a dictionary of features for one cell at one genomic position
        stat_dict = {}
        if (len(rlist) != 0):
            stat_dict['ref_base'] = self.refineBase(rlist[2])
            stat_dict['depth'] = rlist[3]
            stat_dict['refcount'] = self.RefCountString(rlist[4])
            (del_ins_rmvd_str, ins_count, del_count, A_cnt, T_cnt,
             G_cnt, C_cnt, N_cnt) = self.Alt_count(rlist[4])
            stat_dict['ins_count'] = ins_count
            stat_dict['del_count'] = del_count
            stat_dict['A_count'] = A_cnt
            stat_dict['T_count'] = T_cnt
            stat_dict['G_count'] = G_cnt
            stat_dict['C_count'] = C_cnt
            stat_dict['N_count'] = N_cnt
            stat_dict['base_calls'] = self.Create_base_call_string(
                del_ins_rmvd_str, rlist[2])
            stat_dict['base_qual_list'] = self.Get_base_qual_list(rlist[5])
        else:
            # Default empty feature row
            stat_dict['depth'] = 0
            stat_dict['refcount'] = 0
            stat_dict['ref_base'] = ''
            stat_dict['ins_count'] = 0
            stat_dict['del_count'] = 0
            stat_dict['A_count'] = 0
            stat_dict['T_count'] = 0
            stat_dict['G_count'] = 0
            stat_dict['C_count'] = 0
            stat_dict['N_count'] = 0
            stat_dict['base_calls'] = 'NULL'
            stat_dict['base_qual_list'] = []
        return stat_dict

    def calc_strand_bias(self, cell_ftr_pos_list, Alt_count):
        # Compute strand bias using Fisher's exact test
        forward_ref_count = 0
        forward_alt_count = 0
        reverse_ref_count = 0
        reverse_alt_count = 0
        for i in range(len(cell_ftr_pos_list)):
            (fr, fa, rr, ra) = cell_ftr_pos_list[i].Store_Strand_Bias_info()
            forward_ref_count += fr
            forward_alt_count += fa
            reverse_ref_count += rr
            reverse_alt_count += ra
        if ((forward_ref_count == 0) & (reverse_ref_count == 0)):
            return (0.0, 1)
        else:
            cont_table = np.array([[forward_ref_count, reverse_ref_count],
                                   [forward_alt_count, reverse_alt_count]])
            (oddsRatio, pval) = stats.fisher_exact(cont_table)
            return (oddsRatio, pval)

    def calc_prior(self, theta, n_cells, flag):
        # Compute prior probability distribution over number of variant alleles
        prior_variant_number = []
        if flag == 1:
            for i in range(0, 2 * n_cells + 1):
                if ((i == 0) | (i == 2 * n_cells)):
                    lst = [1.0 / i for i in range(1, 2 * n_cells)]
                    prob = 0.5 * (1 - (theta * sum(lst)))
                else:
                    prob = theta / i
                prior_variant_number.append(prob)
        elif flag == 2:
            norm_const = 0
            for i in range(1, 2 * n_cells):
                if (i == n_cells):
                    norm_const += 2 * n_cells
                else:
                    norm_const += float(n_cells) / abs(n_cells - i)
            for i in range(1, 2 * n_cells):
                if (i == n_cells):
                    prob = 2 * n_cells * theta / norm_const
                else:
                    prob = ((n_cells * theta) / abs(n_cells - i)) / norm_const
                prior_variant_number.append(prob)
            sp = sum(prior_variant_number)
            p_l0 = 0.5 * (1 - sp)
            p_l2n = 0.5 * (1 - sp)
            prior_variant_number.insert(0, p_l0)
            prior_variant_number.append(p_l2n)
        elif flag == 3:
            for i in range(0, 2 * n_cells + 1):
                prob = 1. / 2 * n_cells + 1
                prior_variant_number.append(prob)
        return prior_variant_number

    def find_max_prob_ratio(self, matrix, dim):
        # Compute log-likelihood ratio between most likely allele count and zero
        l_0_prob = matrix.denom_prob_matrix[0, dim[1] - 1]
        subtracting_max_prob_val = -743.7469 if l_0_prob == 0 else math.log(l_0_prob)
        (first_term_max_prob_val, max_prob_allele_count) = max(
            (v, i) for i, v in enumerate(matrix.denom_prob_matrix[:, dim[1] - 1]))
        log_first_term_max_prob_val = -743.7469 if first_term_max_prob_val <= 0 else math.log(first_term_max_prob_val)
        max_prob_ratio = log_first_term_max_prob_val - subtracting_max_prob_val
        return (max_prob_ratio, max_prob_allele_count)

    def Get_prior_allele_mat(self, read_smpl_count, alt_smpl_count, cell_no_threshold, total_depth, Alt_freq, pe):
        # Select allele prior based on coverage and allele frequency heuristics
        if ((read_smpl_count > cell_no_threshold - 1) & (alt_smpl_count == 1)):
            prior_allele_mat = allele_prior(0.2)
        elif ((read_smpl_count > cell_no_threshold) & (alt_smpl_count == 2)
              & (total_depth > 30) & (Alt_freq < 0.1)):
            prior_allele_mat = allele_prior(0.1)
        else:
            prior_allele_mat = allele_prior(pe)
        return prior_allele_mat

    def Calc_chr_count(self, barcode):
        # Compute allele count (AC), allele number (AN), and allele frequency (AF)
        AC = 0
        AN = 0
        for c in barcode[1:-1]:
            if c == 'X':
                continue
            else:
                AN += 2
                AC += int(c)
        AF = float(AC) / AN
        return (AC, AF, AN)

    def Calc_Base_Q_Rank_Sum(self, read_supported_cell_list, refBase, altBase):
        # Perform Wilcoxon rank-sum test on base qualities of ref vs alt reads
        ref_read_list = []
        alt_read_list = []
        for cell_ftr_info in read_supported_cell_list:
            for i in range(len(cell_ftr_info.final_bases)):
                if cell_ftr_info.final_bases[i] == refBase:
                    ref_read_list.append(cell_ftr_info.base_qual_int_val_list[i])
                elif cell_ftr_info.final_bases[i] == altBase:
                    alt_read_list.append(cell_ftr_info.base_qual_int_val_list[i])
        return stats.ranksums(alt_read_list, ref_read_list)

    def Calc_qual_depth(self, barcode, All_single_cell_ftrs_list, Qual):
        # Normalize quality by total depth across non-reference genotypes
        depth = 0
        for i in range(len(barcode[1:-1])):
            if barcode[i + 1] == 'X' or barcode[i + 1] == '0':
                continue
            else:
                depth += All_single_cell_ftrs_list[i].depth
        qual_depth = float(Qual) / depth if depth > 0 else Qual
        return qual_depth

    def Get_BAM_RG(self, bam_file):
        # Extract read group ID from BAM header, fallback to filename
        rows = pysam.view("-H", bam_file)
        for r in rows:
            if r.startswith('@RG'):
                r_l = r.split('\t')
                id = r_l[1].split(':')
                return id[1]
        bam_id = bam_file.split('/')[-1].replace('..', '').replace('~', '')
        return bam_id

    def Calc_Per_Smpl_Alt_Ref_Ratio(self, total_ref_depth, Alt_count, read_smpl_count, alt_smpl_count):
        # Compute per-sample normalized alt/ref ratio
        denom = 1 if total_ref_depth == 0 else float(total_ref_depth) / read_smpl_count
        num = float(Alt_count) / alt_smpl_count
        return num / denom

    def Consensus_Filter(self, barcode):
        # Filter variants supported by more than one non-reference genotype
        nu_barcode = barcode[1:-1]
        g_count = 0
        for i in nu_barcode:
            if i in ['1', '2']:
                g_count += 1
        return 1 if g_count > 1 else 0
