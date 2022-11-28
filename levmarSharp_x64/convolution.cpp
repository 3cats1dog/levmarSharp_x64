#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fftw3.h"
#include "convolution.h";

void CplxPlus(fftw_complex res, fftw_complex ab, fftw_complex cd)
{
	res[0] = ab[0] + cd[0];
	res[1] = ab[1] + cd[1];
}
void CplxMinus(fftw_complex res, fftw_complex ab, fftw_complex cd)
{
	res[0] = ab[0] - cd[0];
	res[1] = ab[1] - cd[1];
}
void CplxMltp(fftw_complex res, fftw_complex ab, fftw_complex cd)
{
	res[0] = ab[0] * cd[0] - ab[1] * cd[1];
	res[1] = ab[0] * cd[1] + ab[1] * cd[0];
}
void CplxDiv(fftw_complex res, fftw_complex ab, fftw_complex cd)
{
	if ((cd[0] == 0 && cd[1] == 0) || (ab[0] == 0 && ab[1] == 0))
	{
		res[0] = 0; res[1] = 0;
		return;
	}
	double delta = cd[0] * cd[0] + cd[1] * cd[1];		//cc+dd
	res[0] = (ab[0] * cd[0] + ab[1] * cd[1]) / delta;	//ac+bd
	res[1] = (ab[1] * cd[0] - ab[0] * cd[1]) / delta;	//bc-ad
}


int convolution(double* g, double* f, double* h, int Nf, int Nh)
{

	register unsigned int i;
	fftw_plan p1, p2, p3;

	int N = Nf + Nh;
	double* hFix, *fFix;
	fftw_complex* outh, * outf, * outgC;

	hFix = (double*)fftw_malloc(sizeof(double) * N);
	fFix = (double*)fftw_malloc(sizeof(double) * N);
	outf = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	outh = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	outgC = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

	for (i = 0; i < Nh; i++) { hFix[i] = h[i]; }
	for (i = Nh; i < N; i++) { hFix[i] = 0.0; }

	for (i = 0; i < Nf; i++) { fFix[i] = f[i]; }
	for (i = Nf; i < N; i++) { fFix[i] = 0.0; }


	p1 = fftw_plan_dft_r2c_1d(N, fFix, outf, FFTW_ESTIMATE);
	fftw_execute(p1);

	p2 = fftw_plan_dft_r2c_1d(N, hFix, outh, FFTW_ESTIMATE);
	fftw_execute(p2);

	for (i = 0; i < N; i++) {
		CplxMltp(outgC[i], outf[i], outh[i]);
	}

	p3 = fftw_plan_dft_c2r_1d(N, outgC, g, FFTW_ESTIMATE);
	fftw_execute(p3);


	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);
	fftw_destroy_plan(p3);
	fftw_free(hFix); fftw_free(fFix);
	fftw_free(outf);fftw_free(outh); fftw_free(outgC);

	return 0;
}
int deconvolution(double* g, double* f, double* h, int Nf, int Nh)
{

	register unsigned int i;
	fftw_plan p1, p2, p3;

	int N = Nf+ Nh;
	double* hFix, *fFix;
	fftw_complex* outh, * outf, * outgC;

	hFix = (double*)fftw_malloc(sizeof(double) * N);
	fFix = (double*)fftw_malloc(sizeof(double) * N);
	outf = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	outh = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	outgC = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

	for (i = 0; i < Nh; i++) { hFix[i] = h[i]; }
	for (i = Nh; i < N; i++) { hFix[i] = 0.0; }

	for (i = 0; i < Nf; i++) { fFix[i] = f[i]; }
	for (i = Nf; i < N; i++) { fFix[i] = 0.0; }

	//p1 = fftw_plan_dft_r2c_1d(N, fFix, outf, FFTW_ESTIMATE);
	p1 = fftw_plan_dft_r2c_1d(N, f, outf, FFTW_ESTIMATE);
	fftw_execute(p1);

	p2 = fftw_plan_dft_r2c_1d(N, hFix, outh, FFTW_ESTIMATE);
	fftw_execute(p2);

	for (i = 0; i < N; i++) {
		CplxDiv(outgC[i], outf[i], outh[i]);
	}

	p3 = fftw_plan_dft_c2r_1d(N, outgC, g, FFTW_ESTIMATE);
	fftw_execute(p3);


	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);
	fftw_destroy_plan(p3);
	fftw_free(hFix); fftw_free(fFix);
	fftw_free(outf); fftw_free(outh); fftw_free(outgC);

	return 0;
}

