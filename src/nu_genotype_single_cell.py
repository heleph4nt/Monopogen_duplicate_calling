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


This module defines the probabilistic model for computing
the genotype probability of a single cell, conditional on
the genotypes and evidence from all other cells.
"""

from alleles_prior import allele_prior
from Single_Cell_Ftrs_Pos import Single_Cell_Ftrs_Pos
from genotype_prob_mat import Genotype_Prob_matrix
import sys
import math


class Single_cell_genotype_records:
    """
    This class computes the posterior probability of a given genotype
    (0/0, 0/1, 1/1) for one specific cell, conditioning on:
    - its own read evidence
    - the aggregated evidence from all other cells
    - a prior on the number of variant alleles in the population
    """

    def __init__(
        self,
        current_cell_ftr_list,
        other_cells_ftr_list,
        curr_cell_residual_mat,
        other_cell_list_len,
        n_cells,
        prior_allele_mat,
        prior_variant_allele
    ):
        # Feature object for the current cell
        self.current_cell_ftr_list = current_cell_ftr_list

        # Feature objects for all other cells
        self.other_cells_ftr_list = other_cells_ftr_list

        # Probability matrix summarizing evidence from other cells
        self.residual_mat = curr_cell_residual_mat

        # Number of other cells with read support
        self.other_cell_list_len = other_cell_list_len

        # Total number of cells
        self.n_cells = n_cells

        # Prior allele probability matrix
        self.prior_allele_mat = prior_allele_mat

        # Prior on number of variant alleles
        self.prior_variant_allele = prior_variant_allele

        # Debug containers (not actively used)
        self.p1 = []
        self.p2 = []


    def current_cell_genotype_prob(self, gt_flag):
        """
        Compute the likelihood of the reads in the current cell
        given a candidate genotype (gt_flag ∈ {0,1,2}).
        """
        return self.current_cell_ftr_list.Prob_Reads_Given_Genotype_Genotyping(gt_flag)


    def nCr(self, n, r, factorial_list):
        """
        Compute the binomial coefficient n choose r using precomputed factorials.
        """
        comb_denom = n - r
        number = factorial_list[n] / (
            factorial_list[r] * factorial_list[comb_denom]
        )
        return number


    def find_coeff(self, n_cells, l, j, nCr_matrix):
        """
        Compute the combinatorial coefficient that links:
        - j alternate alleles in the current cell
        - l alternate alleles among all other cells

        This corresponds to:
        (l choose j) * (2n - l choose 2 - j) / (2n choose 2)
        """
        if j > l:
            return 0
        else:
            num = (
                nCr_matrix[l, j] *
                nCr_matrix[2 * n_cells - l, 2 - j]
            )
            denom = nCr_matrix[2 * n_cells, 2]
            return float(num) / denom


    def other_cells_genotype_prob(self, gt_flag, nCr_matrix):
        """
        Compute the probability contribution from all other cells,
        marginalizing over the possible number of alternate alleles
        they may collectively carry.
        """
        if self.other_cell_list_len == 0:
            # If there are no other cells, their contribution is neutral
            return 1
        else:
            prob = 0.0

            # Iterate over possible total alternate allele counts (l)
            for l in range(0, self.residual_mat.dim[0]):
                coeff = self.find_coeff(
                    self.n_cells, l, gt_flag, nCr_matrix
                )

                # Combine:
                # - combinatorial coefficient
                # - probability mass from residual matrix
                # - prior on variant allele counts
                prob += (
                    coeff *
                    self.residual_mat.denom_prob_matrix[
                        l, self.residual_mat.dim[1] - 1
                    ] *
                    self.prior_variant_allele[l + gt_flag]
                )

            return prob


    def find_genotype_prob(self, gt_flag, nCr_matrix):
        """
        Compute the full posterior probability of genotype gt_flag
        for the current cell by multiplying:
        - likelihood from current cell reads
        - marginalized likelihood from other cells
        """
        p1 = self.current_cell_genotype_prob(gt_flag)
        p2 = self.other_cells_genotype_prob(gt_flag, nCr_matrix)

        # Numerical stability safeguards
        if p2 == 0:
            p2 = 1e-322
        if p1 == 0:
            p1 = 1e-322

        return p1 * p2
