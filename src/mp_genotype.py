"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Hamim Zafar and Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

This module handles per-cell genotype inference using multiprocessing.
It computes genotype probabilities for each single cell independently,
conditional on the evidence from all other cells.
"""

from genotype_prob_mat import Genotype_Prob_matrix
from nu_genotype_single_cell import Single_cell_genotype_records
from utils import Utils_Functions
import math
import heapq
import copy
import re
import sys
import numpy as np
from operator import add
import fileinput
from contextlib import closing
from base_q_ascii import Ascii_Table

# Utility functions used throughout the model
U = Utils_Functions()


class MP_single_cell_genotype:
    """
    This class computes genotype likelihoods and genotype calls
    for individual cells, given the variant context and evidence
    from all other cells.

    It is designed to be used with multiprocessing.
    """

    def get_info_string(
        self,
        read_supported_cell_list,
        prior_allele_mat,
        n_cells,
        nCr_matrix,
        prior_variant_number,
        denominator,
        genotype_dict,
        cell_count
    ):
        """
        Compute the genotype call and FORMAT field for one single cell.

        This function:
        - removes the current cell from the population
        - recomputes the probability distribution for the remaining cells
        - computes posterior genotype probabilities for the current cell
        - formats GT, AD, DP, GQ, and PL fields for VCF output
        """

        # Extract the current cell's feature object
        current_cell_ftr_info = read_supported_cell_list[cell_count]

        # Create a list of all other cells except the current one
        (cp_read_supported_cell_list, other_cell_list_len) = U.copy_list_but_one(
            read_supported_cell_list, cell_count
        )

        # Build genotype probability matrix for the remaining cells
        curr_cell_residual_mat = Genotype_Prob_matrix(
            cp_read_supported_cell_list,
            prior_allele_mat,
            other_cell_list_len
        )

        # Fill the probability matrix if there are remaining cells
        if other_cell_list_len > 0:
            curr_cell_residual_mat.denom_prob_matrix = (
                curr_cell_residual_mat.fill_matrix(
                    cp_read_supported_cell_list,
                    n_cells - 1,
                    nCr_matrix
                )
            )

        curr_cell_residual_mat.dim = curr_cell_residual_mat.denom_prob_matrix.shape

        # Create genotype record object for this cell
        cell_genotype_obj = Single_cell_genotype_records(
            current_cell_ftr_info,
            cp_read_supported_cell_list,
            curr_cell_residual_mat,
            other_cell_list_len,
            n_cells,
            prior_allele_mat,
            prior_variant_number
        )

        # Compute posterior genotype probabilities for genotypes 0/0, 0/1, 1/1
        p_g0 = cell_genotype_obj.find_genotype_prob(0, nCr_matrix) / denominator
        p_g1 = cell_genotype_obj.find_genotype_prob(1, nCr_matrix) / denominator
        p_g2 = cell_genotype_obj.find_genotype_prob(2, nCr_matrix) / denominator

        p_list = [p_g0, p_g1, p_g2]

        # Handle degenerate cases where probabilities collapse to zero
        if max(p_list) == 0:
            p_list[0] = current_cell_ftr_info.cell_prob_0
            p_list[1] = current_cell_ftr_info.cell_prob_1
            p_list[2] = current_cell_ftr_info.cell_prob_2

        # Resolve ties using single-cell likelihoods
        if p_list.count(max(p_list)) > 1:
            p_list[0] = current_cell_ftr_info.cell_prob_0
            p_list[1] = current_cell_ftr_info.cell_prob_1
            p_list[2] = current_cell_ftr_info.cell_prob_2

        # Final fallback based on allele frequency
        if max(p_list) == 0:
            if current_cell_ftr_info.Alt_freq < 0.1:
                g_ind = 0
            elif current_cell_ftr_info.Alt_freq > 0.8:
                g_ind = 2
            else:
                g_ind = 1
        else:
            _, g_ind = max((v, i) for i, v in enumerate(p_list))

        # Store genotype index in the cell object
        read_supported_cell_list[cell_count].GT = g_ind

        # Normalize genotype probabilities
        norm_p_list = [0, 0, 0]
        sum_p_list = sum(p_list)
        if sum_p_list == 0:
            norm_p_list[g_ind] = 1
        else:
            for i in range(3):
                norm_p_list[i] = p_list[i] / sum_p_list

        # Compute PL (Phred-scaled genotype likelihoods)
        PL = []
        for p in norm_p_list:
            if p == 0:
                PL.append(3240)
            else:
                PL.append(int(round(-10 * math.log10(p))))

        # Compute GQ (Genotype Quality)
        max_p_g = max(norm_p_list)
        if max_p_g < 1:
            GQ = int(round(-10 * math.log10(1 - max_p_g)))
        else:
            GQ = 3240

        # Map genotype index to GT string
        final_genotype = genotype_dict[g_ind]

        # Encode genotype as barcode character
        cell_barcode = str(g_ind)

        # Build FORMAT field components
        cell_info_list = [
            final_genotype,
            ','.join([
                str(current_cell_ftr_info.refDepth),
                str(current_cell_ftr_info.altcount)
            ]),
            str(current_cell_ftr_info.depth),
            str(GQ),
            ','.join(str(i) for i in PL)
        ]

        cell_info_string = ':'.join(cell_info_list)

        # Return FORMAT string and barcode symbol
        return (cell_info_string, cell_barcode)


    def pre_compute_likelihood(
        self,
        refBase,
        altBase,
        Alt_freq,
        prior_allele_mat,
        max_depth,
        sngle_cell_obj
    ):
        """
        Pre-compute genotype likelihoods for a single cell.
        This function is defined for extensibility but is not
        heavily used in the main Monopogen pipeline.
        """

        cell_prob_0 = sngle_cell_obj.Prob_Reads_Given_Genotype(0, max_depth)
        cell_prob_2 = sngle_cell_obj.Prob_Reads_Given_Genotype(2, max_depth)
        cell_prob_1 = sngle_cell_obj.Prob_Reads_Given_Genotype(1, max_depth)

        return sngle_cell_obj
