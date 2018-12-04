/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef LIB_MATRIX_H_INCLUDED
#define LIB_MATRIX_H_INCLUDED

#include <vector>
#include <string>

/**
 * @brief Multiply two matrix A * B. It uses BLAS SGEMM.
 *
 * @param o_AB = array containing n by l product matrix at exit;
 * @param i_A = input array containing n by m matrix;
 * @param i_B = input array containing m by l matrix;
 * @param p_n, p_m, p_l = dimension parameters of arrays.
 * @param p_transA = true for transposing A.
 * @param p_transA = true for transposing B.
 * @param p_colMajor = true for if matrices should be read by columns.
 *
 * @return  none.
 **/
void productMatrix(
	std::vector<float> &o_AB
,	std::vector<float> const& i_A
,	std::vector<float> const& i_B
,	const unsigned p_n
,	const unsigned p_m
,	const unsigned p_l
,	const bool p_transA
,	const bool p_transB
,	const bool p_colMajor = true
,	unsigned lda = 0
,	unsigned ldb = 0
);

/**
 * @brief Compute a specified number of eigenvectors and eigenvalues of a
 * symmetric matrix.
 *
 * NOTES:
 * - matrices are stored in column-major ordering
 * - columns of input matrices are contiguous in memory
 * - only the upper triangular triangular part of o_mat is used
 * - the upper triangular part of o_mat is destroyed
 * - the output o_U contains the eigenvectors as columns, and is
 *   stored ini column-major ordering (i.e. it returns the eigenvectors
 *   as rows in row-major ordering)
 *
 * @param i_mat: contains input matrix;
 * @param p_n  : size of the matrix;
 * @param p_r  : number of eigenvectors and eigenvalues.
 * @param o_S  : vector with the r eigenvalues
 * @param o_U  : matrix with the r eigenvectors
 *
 * @return none.
 **/
int matrixEigs(
	std::vector<float> &i_mat
,	const unsigned p_n
,	const unsigned p_r
,	std::vector<float> &o_S
,	std::vector<float> &o_U
);

#endif // LIB_MATRIX_H_INCLUDED
