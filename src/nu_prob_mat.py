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


This module builds probability matrices that aggregate read-level
genotype likelihoods across multiple single cells.

These matrices are later used to:
- compute the probability that no variant exists at a site
- marginalize over the number of alternate alleles in the population
"""

import numpy as np
from Single_Cell_Ftrs_Pos import Single_Cell_Ftrs_Pos
from alleles_prior import allele_prior
import math


class Prob_matrix:
    """
    This class constructs a dynamic programming matrix that stores
    the probability of observing the sequencing data given that
    exactly l alternate alleles exist among the cells.

    Rows correspond to total alternate allele counts (l),
    columns correspond to progressively adding cells.
    """

    def __init__(self, sngl_cell_ftr_list, prior_allele_mat, read_supported_n_cells):
        # Number of cells with read support
        self.n_cells = read_supported_n_cells

        # Probability matrix of shape (2*n_cells+1, n_cells)
        self.denom_prob_matrix = np.zeros((2 * self.n_cells + 1, self.n_cells))

        # List of Single_Cell_Ftrs_Pos objects
        self.sngl_cell_ftr_list = sngl_cell_ftr_list

        # Prior allele probability matrix
        self.prior_allele_mat = prior_allele_mat


    def nCr(self, n, r, factorial_list):
        """
        Compute the binomial coefficient (n choose r).
        """
        comb_denom = n - r
        number = factorial_list[n] / (
            factorial_list[r] * factorial_list[comb_denom]
        )
        return number


    def fill_matrix(self, sngl_cell_ftr_list, original_n_cells, max_depth, nCr_matrix, pad):
        """
        Fill the probability matrix using a dynamic programming approach.

        Each entry (l, j) represents the probability of observing the data
        from the first j+1 cells given that l alternate alleles are present.
        """

        # Initialize matrix using the first cell
        self.denom_prob_matrix[0, 0] = float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype(0, max_depth, pad)
        )

        self.denom_prob_matrix[2, 0] = float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype(2, max_depth, pad)
        )

        self.denom_prob_matrix[1, 0] = 2 * float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype(1, max_depth, pad)
        )

        # Iteratively add remaining cells
        for j in range(1, self.n_cells):
            cell_j_prob_0 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype(
                0, max_depth, pad)
            cell_j_prob_2 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype(
                2, max_depth, pad)
            cell_j_prob_1 = 2 * sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype(
                1, max_depth, pad)

            for l in range(0, 2 * self.n_cells + 1):

                # Impossible allele counts are set to zero
                if l > 2 * (j + 1):
                    self.denom_prob_matrix[l, j] = 0
                else:
                    # Transition cases:
                    # - add 0 alt alleles (homozygous reference)
                    # - add 1 alt allele (heterozygous)
                    # - add 2 alt alleles (homozygous alternate)
                    if l == 0:
                        t1 = 0
                        t2 = 0
                        t3 = self.denom_prob_matrix[l, j - 1]
                    elif l == 1:
                        t1 = 0
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]
                    else:
                        t1 = self.denom_prob_matrix[l - 2, j - 1]
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]

                    self.denom_prob_matrix[l, j] = (
                        t1 * cell_j_prob_2 +
                        t2 * cell_j_prob_1 +
                        t3 * cell_j_prob_0
                    )

        # Normalize by combinatorial factor
        for l in range(0, 2 * self.n_cells + 1):
            self.denom_prob_matrix[l, self.n_cells - 1] /= (
                nCr_matrix[2 * self.n_cells, l]
            )

        return self.denom_prob_matrix


    def fill_matrix_50d(self, sngl_cell_ftr_list, original_n_cells, nCr_matrix):
        """
        Alternative matrix fill using a fixed-depth (50 reads) likelihood model.
        """

        self.denom_prob_matrix[0, 0] = float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_50d(0)
        )
        self.denom_prob_matrix[2, 0] = float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_50d(2)
        )
        self.denom_prob_matrix[1, 0] = 2 * float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_50d(1)
        )

        for j in range(1, self.n_cells):
            cell_j_prob_0 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_50d(0)
            cell_j_prob_2 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_50d(2)
            cell_j_prob_1 = 2 * sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_50d(1)

            for l in range(0, 2 * self.n_cells + 1):
                if l > 2 * (j + 1):
                    self.denom_prob_matrix[l, j] = 0
                else:
                    if l == 0:
                        t1 = 0
                        t2 = 0
                        t3 = self.denom_prob_matrix[l, j - 1]
                    elif l == 1:
                        t1 = 0
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]
                    else:
                        t1 = self.denom_prob_matrix[l - 2, j - 1]
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]

                    self.denom_prob_matrix[l, j] = (
                        t1 * cell_j_prob_2 +
                        t2 * cell_j_prob_1 +
                        t3 * cell_j_prob_0
                    )

        for l in range(0, 2 * self.n_cells + 1):
            self.denom_prob_matrix[l, self.n_cells - 1] /= (
                nCr_matrix[2 * self.n_cells, l]
            )

        return self.denom_prob_matrix


    def fill_matrix_stable(self, sngl_cell_ftr_list, original_n_cells, nCr_matrix):
        """
        Numerically stable version of matrix construction.
        Each column is normalized during construction to prevent underflow.
        """

        # Initialize first cell
        self.denom_prob_matrix[0, 0] = sngl_cell_ftr_list[
            0].Prob_Reads_Given_Genotype_Genotyping(0)
        self.denom_prob_matrix[2, 0] = sngl_cell_ftr_list[
            0].Prob_Reads_Given_Genotype_Genotyping(2)
        self.denom_prob_matrix[1, 0] = 2 * sngl_cell_ftr_list[
            0].Prob_Reads_Given_Genotype_Genotyping(1)

        # Normalize first column
        sum_l = (
            self.denom_prob_matrix[0, 0] +
            self.denom_prob_matrix[1, 0] +
            self.denom_prob_matrix[2, 0]
        )
        self.denom_prob_matrix[:, 0] /= sum_l

        # Iteratively add remaining cells with normalization
        for j in range(1, self.n_cells):
            cell_j_prob_0 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(0)
            cell_j_prob_2 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(2)
            cell_j_prob_1 = 2 * sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(1)

            sum_l = 0
            for l in range(0, 2 * self.n_cells + 1):
                if l > 2 * (j + 1):
                    self.denom_prob_matrix[l, j] = 0
                else:
                    if l == 0:
                        t1 = 0
                        t2 = 0
                        t3 = self.denom_prob_matrix[l, j - 1]
                    elif l == 1:
                        t1 = 0
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]
                    else:
                        t1 = self.denom_prob_matrix[l - 2, j - 1]
                        t2 = 2 * self.denom_prob_matrix[l - 1, j - 1]
                        t3 = self.denom_prob_matrix[l, j - 1]

                    self.denom_prob_matrix[l, j] = (
                        t1 * cell_j_prob_2 +
                        t2 * cell_j_prob_1 +
                        t3 * cell_j_prob_0
                    )
                    sum_l += self.denom_prob_matrix[l, j]

            # Normalize column
            for l in range(0, 2 * (j + 1)):
                self.denom_prob_matrix[l, j] /= sum_l

        # Final combinatorial normalization
        for l in range(0, 2 * self.n_cells + 1):
            self.denom_prob_matrix[l, self.n_cells - 1] /= (
                nCr_matrix[2 * self.n_cells, l]
            )

        return self.denom_prob_matrix
