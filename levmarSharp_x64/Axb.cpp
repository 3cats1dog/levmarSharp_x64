/////////////////////////////////////////////////////////////////////////////////
// 
//  Solution of linear systems involved in the Levenberg - Marquardt
//  minimization algorithm
//  Copyright (C) 2004  Manolis Lourakis (lourakis at ics forth gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////////

/******************************************************************************** 
 * LAPACK-based implementations for various linear system solvers. The same core
 * code is used with appropriate #defines to derive single and double precision
 * solver versions, see also Axb_core.c
 ********************************************************************************/
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "levmar.h"
#include "misc.h"

#if !defined(LM_DBL_PREC) && !defined(LM_SNGL_PREC)
#error At least one of LM_DBL_PREC, LM_SNGL_PREC should be defined!
#endif

#ifdef LINSOLVERS_RETAIN_MEMORY
#define __STATIC__ static
#else
#define __STATIC__ // empty
#endif /* LINSOLVERS_RETAIN_MEMORY */


#ifdef LM_DBL_PREC
/* double precision definitions */
#define LM_REAL double
#define LM_CNST(x) (x)
#include <float.h>
#define LM_REAL_EPSILON DBL_EPSILON

int dAx_eq_b_LU_noLapack(LM_REAL* A, LM_REAL* B, LM_REAL* x, int m)
{
	__STATIC__ void* buf = NULL;
	__STATIC__ int buf_sz = 0;

	register int i, j, k;
	int* idx, maxi = -1, idx_sz, a_sz, work_sz, tot_sz;
	LM_REAL* a, * work, max, sum, tmp;

	if (!A)
#ifdef LINSOLVERS_RETAIN_MEMORY
	{
		if (buf) free(buf);
		buf = NULL;
		buf_sz = 0;

		return 1;
	}
#else
		return 1; /* NOP */
#endif /* LINSOLVERS_RETAIN_MEMORY */

  /* calculate required memory size */
	idx_sz = m;
	a_sz = m * m;
	work_sz = m;
	tot_sz = (a_sz + work_sz) * sizeof(LM_REAL) + idx_sz * sizeof(int); /* should be arranged in that order for proper doubles alignment */

#ifdef LINSOLVERS_RETAIN_MEMORY
	if (tot_sz > buf_sz) { /* insufficient memory, allocate a "big" memory chunk at once */
		if (buf) free(buf); /* free previously allocated memory */

		buf_sz = tot_sz;
		buf = (void*)malloc(tot_sz);
		if (!buf) {
			fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_LU) "() failed!\n");
			exit(1);
		}
	}
#else
	buf_sz = tot_sz;
	buf = (void*)malloc(tot_sz);
	if (!buf) {
		fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_LU) "() failed!\n");
		exit(1);
	}
#endif /* LINSOLVERS_RETAIN_MEMORY */

	a = (LM_REAL*)buf;
	work = a + a_sz;
	idx = (int*)(work + work_sz);

	/* avoid destroying A, B by copying them to a, x resp. */
	memcpy(a, A, a_sz * sizeof(LM_REAL));
	memcpy(x, B, m * sizeof(LM_REAL));

	/* compute the LU decomposition of a row permutation of matrix a; the permutation itself is saved in idx[] */
	for (i = 0; i < m; ++i) {
		max = 0.0;
		for (j = 0; j < m; ++j)
			if ((tmp = FABS(a[i * m + j])) > max)
				max = tmp;
		if (max == 0.0) {
			fprintf(stderr, RCAT("Singular matrix A in ", AX_EQ_B_LU) "()!\n");
#ifndef LINSOLVERS_RETAIN_MEMORY
			free(buf);
#endif

			return 0;
		}
		work[i] = LM_CNST(1.0) / max;
	}

	for (j = 0; j < m; ++j) {
		for (i = 0; i < j; ++i) {
			sum = a[i * m + j];
			for (k = 0; k < i; ++k)
				sum -= a[i * m + k] * a[k * m + j];
			a[i * m + j] = sum;
		}
		max = 0.0;
		for (i = j; i < m; ++i) {
			sum = a[i * m + j];
			for (k = 0; k < j; ++k)
				sum -= a[i * m + k] * a[k * m + j];
			a[i * m + j] = sum;
			if ((tmp = work[i] * FABS(sum)) >= max) {
				max = tmp;
				maxi = i;
			}
		}
		if (j != maxi) {
			for (k = 0; k < m; ++k) {
				tmp = a[maxi * m + k];
				a[maxi * m + k] = a[j * m + k];
				a[j * m + k] = tmp;
			}
			work[maxi] = work[j];
		}
		idx[j] = maxi;
		if (a[j * m + j] == 0.0)
			a[j * m + j] = LM_REAL_EPSILON;
		if (j != m - 1) {
			tmp = LM_CNST(1.0) / (a[j * m + j]);
			for (i = j + 1; i < m; ++i)
				a[i * m + j] *= tmp;
		}
	}

	/* The decomposition has now replaced a. Solve the linear system using
	 * forward and back substitution
	 */
	for (i = k = 0; i < m; ++i) {
		j = idx[i];
		sum = x[j];
		x[j] = x[i];
		if (k != 0)
			for (j = k - 1; j < i; ++j)
				sum -= a[i * m + j] * x[j];
		else
			if (sum != 0.0)
				k = i + 1;
		x[i] = sum;
	}

	for (i = m - 1; i >= 0; --i) {
		sum = x[i];
		for (j = i + 1; j < m; ++j)
			sum -= a[i * m + j] * x[j];
		x[i] = sum / a[i * m + i];
	}

#ifndef LINSOLVERS_RETAIN_MEMORY
	free(buf);
#endif

	return 1;
}


int dAux_conv_noLapack(LM_REAL* G, LM_REAL* X, LM_REAL* h, int n, int m)
{
	__STATIC__ void* buf = NULL;
	__STATIC__ int buf_sz = 0;

	register int i, j, js;
	int diff_sz, tot_sz;
	LM_REAL * diff, max, sum, tmp;

	if (!G)
#ifdef LINSOLVERS_RETAIN_MEMORY
	{
		if (buf) free(buf);
		buf = NULL;
		buf_sz = 0;

		return 1;
	}
#else
		return 1; /* NOP */
#endif /* LINSOLVERS_RETAIN_MEMORY */

  /* calculate required memory size */
	diff_sz = n;
	tot_sz = (diff_sz) * sizeof(LM_REAL);  /* should be arranged in that order for proper doubles alignment */

#ifdef LINSOLVERS_RETAIN_MEMORY
	if (tot_sz > buf_sz) { /* insufficient memory, allocate a "big" memory chunk at once */
		if (buf) free(buf); /* free previously allocated memory */

		buf_sz = tot_sz;
		buf = (void*)malloc(tot_sz);
		if (!buf) {
			fprintf(stderr, RCAT("memory allocation in dAux_conv_noLapack failed!\n"));
			exit(1);
		}
	}
#else
	buf_sz = tot_sz;
	buf = (void*)malloc(tot_sz);
	if (!buf) {
		fprintf(stderr, RCAT("memory allocation in dAux_conv_noLapack failed!\n"));
		exit(1);
	}
#endif /* LINSOLVERS_RETAIN_MEMORY */

	diff = (LM_REAL*)buf;

	diff[0] = 0;
	for (i = 1; i < diff_sz; i++)
	{
		diff[i] = (LM_REAL)(X[i] - X[i - 1]);
	}

	for (i = 0; i < n; i++)
	{
		G[i] = 0;
		js = 0;
		if (i > m)
		{
			js = (i - m + 1);
			G[i] += X[js] - X[0];
		}
		for (j = js; j <= i; j++)
		{
			G[i] += diff[j] * h[i - j];
		}
	}

#ifndef LINSOLVERS_RETAIN_MEMORY
	free(buf);
#endif

	return 1;
}

#undef LM_REAL
#undef LM_CNST
#undef LM_REAL_EPSILON
#endif /* LM_DBL_PREC */

#ifdef LM_SNGL_PREC
/* single precision (float) definitions */
#define LM_REAL float
#define __SUBCNST(x) x##F
#define LM_CNST(x) __SUBCNST(x) // force substitution
#define LM_REAL_EPSILON FLT_EPSILON

int sAx_eq_b_LU_noLapack(LM_REAL* A, LM_REAL* B, LM_REAL* x, int m)
{
	__STATIC__ void* buf = NULL;
	__STATIC__ int buf_sz = 0;

	register int i, j, k;
	int* idx, maxi = -1, idx_sz, a_sz, work_sz, tot_sz;
	LM_REAL* a, * work, max, sum, tmp;

	if (!A)
#ifdef LINSOLVERS_RETAIN_MEMORY
	{
		if (buf) free(buf);
		buf = NULL;
		buf_sz = 0;

		return 1;
	}
#else
		return 1; /* NOP */
#endif /* LINSOLVERS_RETAIN_MEMORY */

  /* calculate required memory size */
	idx_sz = m;
	a_sz = m * m;
	work_sz = m;
	tot_sz = (a_sz + work_sz) * sizeof(LM_REAL) + idx_sz * sizeof(int); /* should be arranged in that order for proper doubles alignment */

#ifdef LINSOLVERS_RETAIN_MEMORY
	if (tot_sz > buf_sz) { /* insufficient memory, allocate a "big" memory chunk at once */
		if (buf) free(buf); /* free previously allocated memory */

		buf_sz = tot_sz;
		buf = (void*)malloc(tot_sz);
		if (!buf) {
			fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_LU) "() failed!\n");
			exit(1);
		}
	}
#else
	buf_sz = tot_sz;
	buf = (void*)malloc(tot_sz);
	if (!buf) {
		fprintf(stderr, RCAT("memory allocation in ", AX_EQ_B_LU) "() failed!\n");
		exit(1);
	}
#endif /* LINSOLVERS_RETAIN_MEMORY */

	a = (LM_REAL*)buf;
	work = a + a_sz;
	idx = (int*)(work + work_sz);

	/* avoid destroying A, B by copying them to a, x resp. */
	memcpy(a, A, a_sz * sizeof(LM_REAL));
	memcpy(x, B, m * sizeof(LM_REAL));

	/* compute the LU decomposition of a row permutation of matrix a; the permutation itself is saved in idx[] */
	for (i = 0; i < m; ++i) {
		max = 0.0;
		for (j = 0; j < m; ++j)
			if ((tmp = FABS(a[i * m + j])) > max)
				max = tmp;
		if (max == 0.0) {
			fprintf(stderr, RCAT("Singular matrix A in ", AX_EQ_B_LU) "()!\n");
#ifndef LINSOLVERS_RETAIN_MEMORY
			free(buf);
#endif

			return 0;
		}
		work[i] = LM_CNST(1.0) / max;
	}

	for (j = 0; j < m; ++j) {
		for (i = 0; i < j; ++i) {
			sum = a[i * m + j];
			for (k = 0; k < i; ++k)
				sum -= a[i * m + k] * a[k * m + j];
			a[i * m + j] = sum;
		}
		max = 0.0;
		for (i = j; i < m; ++i) {
			sum = a[i * m + j];
			for (k = 0; k < j; ++k)
				sum -= a[i * m + k] * a[k * m + j];
			a[i * m + j] = sum;
			if ((tmp = work[i] * FABS(sum)) >= max) {
				max = tmp;
				maxi = i;
			}
		}
		if (j != maxi) {
			for (k = 0; k < m; ++k) {
				tmp = a[maxi * m + k];
				a[maxi * m + k] = a[j * m + k];
				a[j * m + k] = tmp;
			}
			work[maxi] = work[j];
		}
		idx[j] = maxi;
		if (a[j * m + j] == 0.0)
			a[j * m + j] = LM_REAL_EPSILON;
		if (j != m - 1) {
			tmp = LM_CNST(1.0) / (a[j * m + j]);
			for (i = j + 1; i < m; ++i)
				a[i * m + j] *= tmp;
		}
	}

	/* The decomposition has now replaced a. Solve the linear system using
	 * forward and back substitution
	 */
	for (i = k = 0; i < m; ++i) {
		j = idx[i];
		sum = x[j];
		x[j] = x[i];
		if (k != 0)
			for (j = k - 1; j < i; ++j)
				sum -= a[i * m + j] * x[j];
		else
			if (sum != 0.0)
				k = i + 1;
		x[i] = sum;
	}

	for (i = m - 1; i >= 0; --i) {
		sum = x[i];
		for (j = i + 1; j < m; ++j)
			sum -= a[i * m + j] * x[j];
		x[i] = sum / a[i * m + i];
	}

#ifndef LINSOLVERS_RETAIN_MEMORY
	free(buf);
#endif

	return 1;
}


#undef LM_REAL
#undef __SUBCNST
#undef LM_CNST
#undef LM_REAL_EPSILON
#endif /* LM_SNGL_PREC */


#undef __STATIC__