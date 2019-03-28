#include "diag.h"

extern "C" {

void dsaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol, void* resid, 
		int* ncv, void* v, int* ldv, int* iparam, int* ipntr, void* workd,
	   	void* workl, int* lworkl, int* info);

void dseupd_(int* rvec, char* A, int* select, void* d, void* z, int* ldz, void* sigma,
		char* bmat, int* n, char* which, int* nev, double* tol, void* resid, int* ncv, 
		void* v, int* ldv, int* iparam, int* ipntr, void* workd, void* workl, int* lworkl, int* ierr);
}


void GenMatProd::perform_op(const MatType* xIn, MatType* yOut)
{
	char uplo = 'U';
	MKL_INT n = baseHam_.end - baseHam_.start;
	mkl_dcsrsymv(&uplo, &n , baseHam_.mat , baseHam_.ia, baseHam_.ja , xIn , yOut);
}


void calcEValues(const SparseMat& baseHam, GenMatProd& op, MatType* evalues, MatType* evecs, int neigs)
{
	Count maxn = baseHam.end - baseHam.start; 	
	int maxncv = 4 * neigs + 1;
	int ldv = maxn;
	int iparam[11] = {0};
	int ipntr[14];
	int select[maxncv];
	MatType* d = new MatType[maxncv];
	MatType* Z = new MatType[ldv * maxncv];
	MatType* v = new MatType[ldv * maxncv];
	MatType* workd = new MatType[3 * maxn];
	MatType* workev = new MatType[2 * maxncv];
	MatType* resid = new MatType[maxn];
	MatType* workl = new MatType[3 * maxncv * maxncv + 5 * maxncv];
	char bmat[2] = "I";
	char which[3] = "SA";
	int ido, n, nx, nev, ncv, lworkl, info, ierr, ishfts, maxitr, mode1;
	MatType sigma = 0.0;
	double tol;
	int rvec;

	nx = maxn;
	n = nx;
	nev = neigs;
	ncv = 3 * neigs + 1;
	lworkl  = 3 * ncv * ncv + 5 * ncv;
	tol    = 1E-6;
	ido    = 0;
	info   = 0;
	
	ishfts = 1;
	maxitr = 500;
	mode1 = 1;
	iparam[0] = ishfts;
	iparam[2] = maxitr;
	iparam[3] = 1;
	iparam[6] = mode1;

	while (true) {
		dsaupd_(&ido, bmat, &n, which, &nev, &tol, (void*) resid, &ncv,
				(void*) v, &ldv, iparam, ipntr, (void*) workd, (void*) workl, &lworkl, &info);


		if (ido == -1 || ido == 1) {
			op.perform_op(&workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
		}
		else {
			break;
		}
	}

	if (info >= 0) {
		rvec = 1;
		char A = 'A';

		dseupd_(&rvec, &A, select, (void*) d, (void*) Z, &ldv, &sigma, bmat, &n, which, &nev, 
				&tol, (void*) resid, &ncv, (void*) v, &ldv, iparam, ipntr, (void*) workd, (void*) workl, &lworkl, &ierr);

		for (int i = neigs - 1; i >= 0; --i) {
			evalues[neigs - i - 1] = d[i];
			for (int j = 0; j < n; ++j) {
				evecs[(neigs - i - 1) * n + j] = Z[i * n + j];
			}
		}
	}

	delete[] d;
	delete[] Z;
	delete[] v;
	delete[] workd;
	delete[] workev;
	delete[] resid;
	delete[] workl;
}

