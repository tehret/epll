/*
 * Modified work Copyright (c) 2014, Pablo Arias <pariasm@gmail.com>
 * Original work Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file LibMatrix.cpp
 * @brief Tools for matrix manipulation, based on ccmath functions
 *        by Daniel A. Atkinson.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "LibMatrix.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <cblas.h>
#include <lapacke.h>

using namespace std;

/**
 * @brief Multiply two matrix A * B. It uses BLAS SGEMM.
 *
 * @param o_AB = array containing nxl product matrix at exit;
 * @param i_A = input array containing nxm (or mxn if p_transA) matrix;
 * @param i_B = input array containing mxl (or lxm if p_transB) matrix;
 * @param p_n, p_l, p_m = dimension parameters of arrays.
 * @param p_transA = true for transposing A.
 * @param p_transA = true for transposing B.
 * @param p_colMajor = true for if matrices should be read by columns.
 *
 * @return  none.
 **/
void productMatrix(
	vector<float> &o_AB
,	vector<float> const& i_A
,	vector<float> const& i_B
,	const unsigned p_n
,	const unsigned p_l
,	const unsigned p_m
,	const bool p_transA
,	const bool p_transB
,	const bool p_colMajor
){
	unsigned lda = p_colMajor ? (p_transA ? p_m : p_n) : (p_transA ? p_n : p_m);
	unsigned ldb = p_colMajor ? (p_transB ? p_l : p_m) : (p_transB ? p_m : p_l);

	cblas_sgemm(p_colMajor ? CblasColMajor : CblasRowMajor,  // matrix storage mode
	            p_transA   ? CblasTrans    : CblasNoTrans,   // op(A)
	            p_transB   ? CblasTrans    : CblasNoTrans,   // op(B)
	            p_n,                                         // rows(op(A)) [= rows(AB)   ]
	            p_l,                                         // cols(op(B)) [= cols(AB)   ]
	            p_m,                                         // cols(op(A)) [= rows(op(B))]
	            1.f,                                         // alpha
	            i_A.data(), lda,                             // A, lda
	            i_B.data(), ldb,                             // B, ldb
	            0.f,                                         // beta
	            o_AB.data(), p_colMajor ? p_n : p_l          // AB, ldab
	            );
}

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
	vector<float> &i_mat
,	const unsigned p_n
,	const unsigned p_r
,	vector<float> &o_S
,	vector<float> &o_U
){
	// set parameters for LAPACKE SSYEV function
	// SSYEVX: Single SYmmetric EigenValues and eigenvectors eXpert
	lapack_int m;          //< total values of eigenvalues found

	o_S.resize(p_n);       //< array of dimension n. The first m entries
	                       //< contain the eigenvalues found in ascending order.

	o_U.resize(p_n*p_r);   //< the first m columns contain the output eigenvectors

	lapack_int lda = p_n;  //< leading dimension for input matrix
	lapack_int ldu = p_n;  //< leading dimension for output eigenvectors matrix

	lapack_int ifail[p_n]; //< if jobz == 'V' and info > 0, then ifail contains
	                       //< the indices of the eigenvectors that didn't converge
	
	lapack_int info;       //< info =  0 : successful exit
	                       //< info = -i : ith argument is wrong
	                       //< info =  i : i eigenvectors failed to converge

	info = LAPACKE_ssyevx(LAPACK_COL_MAJOR,
			'V',                // compute eigenvalues and eigenvectors
			'I',                // range 'I': eigenvals/vecs between IL-th and IU-th
			'U',                // use upper triangular part of A
			p_n,                // order of matrix A
			i_mat.data(),       // matrix A
			lda,                // stride of matrix
			-1, -1,             // not used (used only when range is 'V'
			p_n - p_r + 1, p_n, // IL and IU indices for eigenvals/vecs range
			0,                  // abstol for stopping criterion
			&m,                 // total values of eigenvectors found
			o_S.data(),         // eigenvalues output
			o_U.data(),         // eigenvectors matrix output
			ldu,                // eigenvectors matrix stride
			ifail               // eigenvectors that did not converge
			);

	return(info);
}
