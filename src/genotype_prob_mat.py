"""
The MIT License
...
"""

import numpy as np
from Single_Cell_Ftrs_Pos import Single_Cell_Ftrs_Pos
from alleles_prior import allele_prior
import math


class Genotype_Prob_matrix:
    # This class builds a probability matrix used to compute the
    # denominator of genotype likelihoods across multiple single cells.
    #
    # Conceptually, it aggregates per-cell genotype likelihoods
    # into a joint distribution over the total number of alternative
    # alleles observed across all cells.
    #
    # This formulation originates from Monovar and is reused by Monopogen
    # to perform probabilistic genotyping under sparse single-cell coverage.

    def __init__(self, sngl_cell_ftr_list, prior_allele_mat, read_supported_n_cells_others):
        # Number of cells contributing reads at this locus
        self.n_cells = read_supported_n_cells_others

        # Probability matrix:
        # rows   = total number of alternative alleles (0 … 2*n_cells)
        # columns = incremental inclusion of cells (0 … n_cells-1)
        self.denom_prob_matrix = np.zeros((2 * self.n_cells + 1, self.n_cells))

        # List of per-cell feature objects (Single_Cell_Ftrs_Pos)
        self.sngl_cell_ftr_list = sngl_cell_ftr_list

        # Prior allele probability matrix (population-level prior)
        self.prior_allele_mat = prior_allele_mat


    def nCr(self, n, r, factorial_list):
        # Compute binomial coefficient "n choose r"
        # Used for normalization over genotype combinations
        comb_denom = n - r
        number = factorial_list[n] / \
            (factorial_list[r] * factorial_list[comb_denom])
        return number


    def fill_matrix(self, sngl_cell_ftr_list, original_n_cells, nCr_matrix):
        # Build the denominator probability matrix using a dynamic programming approach.
        #
        # Each column j incorporates one additional cell.
        # Each row l represents the total number of alternative alleles
        # accumulated across the processed cells.

        # Initialize matrix for the first cell
        self.denom_prob_matrix[0, 0] = float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(0))

        self.denom_prob_matrix[2, 0] = float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(2))

        # Heterozygous genotype has two equivalent allele configurations
        self.denom_prob_matrix[1, 0] = 2 * float(
            sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(1))

        # Iteratively incorporate remaining cells
        for j in range(1, self.n_cells):
            cell_j_prob_0 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(0)
            cell_j_prob_2 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(2)
            cell_j_prob_1 = 2 * sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(1)

            for l in range(0, 2 * self.n_cells + 1):

                # Impossible allele counts
                if l > 2 * (j + 1):
                    self.denom_prob_matrix[l, j] = 0
                else:
                    # Transition cases:
                    # t1: previous state contributes two alt alleles
                    # t2: previous state contributes one alt allele
                    # t3: previous state contributes zero alt alleles
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

                    # Combine probabilities from all genotype possibilities
                    self.denom_prob_matrix[l, j] = (
                        t1 * cell_j_prob_2 +
                        t2 * cell_j_prob_1 +
                        t3 * cell_j_prob_0
                    )

        # Normalize final column using combinatorial counts
        for l in range(0, 2 * self.n_cells + 1):
            self.denom_prob_matrix[l, self.n_cells - 1] /= nCr_matrix[2 * self.n_cells, l]

        return self.denom_prob_matrix


    def fill_matrix_stable(self, sngl_cell_ftr_list, original_n_cells, nCr_matrix):
        # Numerically stable version of fill_matrix().
        #
        # Probabilities are normalized at each step to avoid underflow,
        # which is critical when dealing with many cells or low read depth.

        # Initialize first cell
        self.denom_prob_matrix[0, 0] = sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(0)
        self.denom_prob_matrix[2, 0] = sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(2)
        self.denom_prob_matrix[1, 0] = 2 * sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(1)

        # Normalize initial column
        sum_l = (
            self.denom_prob_matrix[0, 0] +
            self.denom_prob_matrix[1, 0] +
            self.denom_prob_matrix[2, 0]
        )
        self.denom_prob_matrix[:, 0] /= sum_l

        # Iteratively process remaining cells
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

            # Normalize column j
            for l in range(0, 2 * (j + 1)):
                self.denom_prob_matrix[l, j] /= sum_l

        # Final combinatorial normalization
        for l in range(0, 2 * self.n_cells + 1):
            self.denom_prob_matrix[l, self.n_cells - 1] /= nCr_matrix[2 * self.n_cells, l]

        return self.denom_prob_matrix
