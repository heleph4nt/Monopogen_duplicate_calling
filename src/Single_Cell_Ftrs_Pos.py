"""
This module defines the Single_Cell_Ftrs_Pos class.

It represents all sequencing-derived features for ONE CELL
at ONE GENOMIC POSITION, and provides methods to compute
likelihoods P(reads | genotype).

This is the lowest-level probabilistic unit in MonoVar / Monopogen.
"""

from alleles_prior import allele_prior
from utils import Utils_Functions
import math
import copy

U = Utils_Functions()


class Single_Cell_Ftrs_Pos:
    """
    This class stores and processes pileup-derived features
    for a single cell at a single genomic position.

    It handles:
    - parsing pileup strings
    - removing indels
    - extracting base qualities
    - computing genotype likelihoods under sequencing error,
      allelic dropout (ADO), and base quality uncertainty
    """

    # Constructor: initialize features from pileup fields
    def __init__(self, refBase, current_pos_info_list):
        # If depth is zero, initialize empty cell
        if (int(current_pos_info_list[0]) == 0):
            self.depth = 0
            self.refDepth = 0
        else:
            self.refBase = refBase
            self.depth = int(current_pos_info_list[0])
            self.primary_bases = current_pos_info_list[1]
            self.base_q = current_pos_info_list[2]

            # Count reference-supporting reads on forward and reverse strands
            (self.forward_ref_count,
             self.reverse_ref_count,
             self.refDepth) = U.RefCountString(self.primary_bases)


    # Remove insertions and deletions from pileup base string
    def Get_Ins_Del_rmvd_bases(self):
        if ((self.primary_bases.count('+') + self.primary_bases.count('-')) == 0):
            # No indels present
            self.ins_count = 0
            self.del_count = 0
            self.ins_list = []
            self.del_list = []
            self.ins_del_rmvd_bases = self.primary_bases
        else:
            # Remove insertions
            cp_primary_bases = copy.deepcopy(self.primary_bases)
            (self.ins_list, ins_rmvd_bases) = U.find_indel(
                cp_primary_bases, '\+[0-9]+[ACGTNacgtn]+')
            self.ins_count = len(self.ins_list)

            # Remove deletions
            (self.del_list, self.ins_del_rmvd_bases) = U.find_indel(
                ins_rmvd_bases, '-[0-9]+[ACGTNacgtn]+')
            self.del_count = len(self.del_list)
        return 0


    # Convert ASCII base quality characters into error probabilities
    def Get_Base_Qual_Vals(self):
        (self.base_qual_val_list,
         self.base_qual_int_val_list) = U.Get_base_qual_list(self.base_q)
        return 0


    # Generate cleaned base-call string after removing indels and read markers
    def Get_Base_Calls(self, ref):
        (self.start_read_counts,
         self.end_read_counts,
         self.start_end_ins_del_rmvd_bases) = U.Count_Start_and_End(
            self.ins_del_rmvd_bases)

        self.final_bases = U.Create_base_call_string(
            self.start_end_ins_del_rmvd_bases, ref)
        return 0


    # Full preprocessing pipeline for bases and qualities
    def Get_base_call_string_nd_quals(self, ref):
        self.Get_Ins_Del_rmvd_bases()
        self.Get_Base_Qual_Vals()
        self.Get_Base_Calls(ref)
        return 0


    # Count occurrences of each nucleotide among cleaned bases
    def Get_Alt_Allele_Count(self):
        self.A_cnt = self.start_end_ins_del_rmvd_bases.count('A') + \
                     self.start_end_ins_del_rmvd_bases.count('a')
        self.C_cnt = self.start_end_ins_del_rmvd_bases.count('C') + \
                     self.start_end_ins_del_rmvd_bases.count('c')
        self.G_cnt = self.start_end_ins_del_rmvd_bases.count('G') + \
                     self.start_end_ins_del_rmvd_bases.count('g')
        self.T_cnt = self.start_end_ins_del_rmvd_bases.count('T') + \
                     self.start_end_ins_del_rmvd_bases.count('t')
        return 0


    # Store index of the cell in the population
    def Set_Cell_Index(self, index):
        self.cell_index = index
        return 0


    # Store information shared across all cells at this site
    def Store_Addl_Info(self, refBase, altBase, Alt_freq, prior_allele_mat):
        self.altBase = altBase
        self.Alt_freq = Alt_freq
        self.prior_allele_mat = prior_allele_mat
        self.refBase = refBase
        return 0


    # Extract strand-specific counts for strand-bias testing
    def Store_Strand_Bias_info(self):
        self.lowercase_alt_base = self.altBase.lower()
        self.forward_alt_count = self.start_end_ins_del_rmvd_bases.count(self.altBase)
        self.reverse_alt_count = self.start_end_ins_del_rmvd_bases.count(self.lowercase_alt_base)
        self.altcount = self.forward_alt_count + self.reverse_alt_count
        return (self.forward_ref_count, self.forward_alt_count,
                self.reverse_ref_count, self.reverse_alt_count)


    # Normalize heterozygous allele ordering (e.g. GA → AG)
    def refineG(self, g):
        if g == 'CA': g = 'AC'
        elif g == 'GA': g = 'AG'
        elif g == 'TA': g = 'AT'
        elif g == 'GC': g = 'CG'
        elif g == 'TC': g = 'CT'
        elif g == 'TG': g = 'GT'
        return g


    # Compute likelihood P(reads | genotype) for a given genotype
    def Calc_Prob_gt(self, gt, max_depth):
        val = 1.0
        ub = min(len(self.base_qual_val_list),
                 len(self.final_bases), max_depth)
        for i in range(ub):
            curr_base = self.final_bases[i]
            curr_base_key = (gt, curr_base)
            curr_err = self.base_qual_val_list[i]
            prob_i = self.prior_allele_mat.getValue(curr_base_key)

            # Sequencing error model
            prob = curr_err * (1 - prob_i) / 3 + (1 - curr_err) * prob_i
            val *= prob
        return val
