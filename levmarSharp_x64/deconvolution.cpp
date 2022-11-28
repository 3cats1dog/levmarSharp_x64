#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
//
//double PI;
//typedef _Dcomplex cplx;
//
//void _fft(cplx buf[], cplx out[], int n, int step)
//{
//	if (step < n) {
//		_fft(out, buf, n, step * 2);
//		_fft(out + step, buf + step, n, step * 2);
//
//		for (int i = 0; i < n; i += 2 * step) {
//			cplx t = cexp(-I * PI * i / n) * out[i + step];
//			buf[i / 2] = out[i] + t;
//			buf[(i + n) / 2] = out[i] - t;
//		}
//	}
//}
//
//void fft(cplx buf[], int n)
//{
//	cplx out[n];
//	for (int i = 0; i < n; i++) out[i] = buf[i];
//	_fft(buf, out, n, 1);
//}
//
///* pad array length to power of two */
//cplx* pad_two(double g[], int len, int* ns)
//{
//	int n = 1;
//	if (*ns) n = *ns;
//	else while (n < len) n *= 2;
//
//	cplx* buf =(cplx *)calloc(sizeof(cplx), n);
//	for (int i = 0; i < len; i++) buf[i]._Val[0] = g[i];
//	*ns = n;
//	return buf;
//}
//
//void deconv(double g[], int lg, double f[], int lf, double out[]) {
//	int ns = 0;
//	cplx* g2 = pad_two(g, lg, &ns);
//	cplx* f2 = pad_two(f, lf, &ns);
//
//	fft(g2, ns);
//	fft(f2, ns);
//
//	cplx h[ns];
//
//	for (int i = 0; i < ns; i++) h[i] = g2[i] / f2[i];
//	fft(h, ns);
//
//	for (int i = 0; i >= lf - lg; i--)
//		out[-i] = h[(i + ns) % ns] / 32;
//	free(g2);
//	free(f2);
//}