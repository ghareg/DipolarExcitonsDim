#include "opermask.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>


int constrDif(const Occup* qstates, Count lftNum, Count rgtNum, Index* lft, Index* rgt, const Param& pm)
{
	if (lftNum == rgtNum) {
		return 0;
	}
	const Occup* lftOccup = &qstates[lftNum * pm.qSize];
	const Occup* rgtOccup = &qstates[rgtNum * pm.qSize];
	int lftCount = 0;
	int rgtCount = 0;
	int dif = 0;
	Index ind = 0;
	while(ind < pm.qSize) {
		if (lftOccup[ind] == rgtOccup[ind]) {
			++ind;
		}
		else {
			if (lftOccup[ind] > rgtOccup[ind]) {
				dif = lftOccup[ind] - rgtOccup[ind];
				if (lftCount + dif > 2) {
					return -1;
				}
				lft[lftCount++] = ind;
				if (dif > 1) {
					lft[lftCount++] = ind;
				}
			}
			else {
				dif = rgtOccup[ind] - lftOccup[ind];	
				if (rgtCount + dif > 2) {
					return -1;
				}
				rgt[rgtCount++] = ind;
				if (dif > 1) {
					rgt[rgtCount++] = ind;
				}
			}
			++ind;
		}
	}

	return lftCount;
}

double MatElement(Count lftNum, Count rgtNum, const State* qnummap, const Occup* qstates, const Param& pm)
{
	Index lft[2];
	Index rgt[2];
	int diffNum = constrDif(qstates, lftNum, rgtNum, lft, rgt, pm);
	if (diffNum == -1) {
	   return 0.0;
	}

	Occup const* lftOccup = &qstates[lftNum * pm.qSize];
	Occup const* rgtOccup = &qstates[rgtNum * pm.qSize];

	double OneMatElem = 0.0;
	if (diffNum == 0) {
		for (int i = 0; i < pm.qSize; ++i) {
			if (lftOccup[i] != 0) {
				OneMatElem += lftOccup[i] * OneBodyElement(qnummap[i], pm);
			}
		}
	}
	
	double TwoMatElem = 0.0;
	if (diffNum == 0) {
		for (int i = 0; i < pm.qSize; ++i) {
			if (lftOccup[i] > 0) {
				for (int j = i + 1; j < pm.qSize; ++j) {
					if (lftOccup[j] > 0) {
						TwoMatElem += lftOccup[i] * lftOccup[j] * (TwoBodyElement(qnummap[i], qnummap[j], qnummap[j], qnummap[i], pm) + TwoBodyElement(qnummap[i], qnummap[j], qnummap[i], qnummap[j], pm));
					}
				}
				if (lftOccup[i] > 1) {
					TwoMatElem += 0.5 * lftOccup[i] * (lftOccup[i] - 1) * TwoBodyElement(qnummap[i], qnummap[i], qnummap[i], qnummap[i], pm);
				}
			}
		}
	}
	else if (diffNum == 1) {
		for (int i = 0; i < pm.qSize; ++i) {
			if (i == lft[0]) {
				if (lftOccup[i] > 1 && rgtOccup[i] > 0) {
					TwoMatElem += 0.5 * sqrt(lftOccup[i] * (lftOccup[i] - 1) * rgtOccup[i] * rgtOccup[rgt[0]]) * (TwoBodyElement(qnummap[i], qnummap[i], qnummap[i], qnummap[rgt[0]], pm) + TwoBodyElement(qnummap[i], qnummap[i], qnummap[rgt[0]], qnummap[i], pm));
				}
			}
			else if (i == rgt[0]) {
				if (lftOccup[i] > 0 && rgtOccup[i] > 1) {
					TwoMatElem += 0.5 * sqrt(lftOccup[i] * lftOccup[lft[0]] * rgtOccup[i] * (rgtOccup[i] - 1)) * (TwoBodyElement(qnummap[i], qnummap[lft[0]], qnummap[i], qnummap[i], pm) + TwoBodyElement(qnummap[lft[0]], qnummap[i], qnummap[i], qnummap[i], pm));
				}
			}
			else {
				if (lftOccup[i] > 0 && rgtOccup[i] > 0) {
					TwoMatElem += sqrt(lftOccup[i] * lftOccup[lft[0]] * rgtOccup[i] * rgtOccup[rgt[0]]) * (TwoBodyElement(qnummap[lft[0]], qnummap[i], qnummap[i], qnummap[rgt[0]], pm) + TwoBodyElement(qnummap[lft[0]], qnummap[i], qnummap[rgt[0]], qnummap[i], pm));
				}
			}
		}
	}
	else {
		if (lft[0] == lft[1] && rgt[0] == rgt[1]) {
			TwoMatElem += 0.5 * sqrt(lftOccup[lft[0]] * (lftOccup[lft[0]] - 1) * rgtOccup[rgt[0]] * (rgtOccup[rgt[0]] - 1)) * TwoBodyElement(qnummap[lft[0]], qnummap[lft[0]], qnummap[rgt[0]], qnummap[rgt[0]], pm);
		}
		else if (lft[0] == lft[1]) {
			TwoMatElem += 0.5 * sqrt(lftOccup[lft[0]] * (lftOccup[lft[0]] - 1) * rgtOccup[rgt[0]] * rgtOccup[rgt[1]]) * (TwoBodyElement(qnummap[lft[0]], qnummap[lft[0]], qnummap[rgt[0]], qnummap[rgt[1]], pm) + TwoBodyElement(qnummap[lft[0]], qnummap[lft[0]], qnummap[rgt[1]], qnummap[rgt[0]], pm));
		}
		else if (rgt[0] == rgt[1]) {
			TwoMatElem += 0.5 * sqrt(lftOccup[lft[0]] * lftOccup[lft[1]] * rgtOccup[rgt[0]] * (rgtOccup[rgt[0]] - 1)) * (TwoBodyElement(qnummap[lft[0]], qnummap[lft[1]], qnummap[rgt[0]], qnummap[rgt[0]], pm) + TwoBodyElement(qnummap[lft[1]], qnummap[lft[0]], qnummap[rgt[0]], qnummap[rgt[0]], pm));
		}
		else {
			TwoMatElem += sqrt(lftOccup[lft[0]] * lftOccup[lft[1]] * rgtOccup[rgt[0]] * rgtOccup[rgt[1]]) * (TwoBodyElement(qnummap[lft[0]], qnummap[lft[1]], qnummap[rgt[0]], qnummap[rgt[1]], pm) + TwoBodyElement(qnummap[lft[0]], qnummap[lft[1]], qnummap[rgt[1]], qnummap[rgt[0]], pm));
		}
	}

	return OneMatElem + TwoMatElem;
}

double OneBodyElement(const State& st, const Param& pm)
{
	return pm.tHop * cos(st.dir * dk + k0);
}

double TwoBodyElement(const State& lft1, const State& lft2, const State& rgt1, const State& rgt2, const Param& pm)
{
	double result = 0.0;
	if (lft1.layer == rgt2.layer || lft2.layer == rgt1.layer) {
	   if (lft1.dir + lft2.dir == rgt1.dir + rgt2.dir) {
		   if (lft1.layer == lft2.layer) {
			   if(!pm.quasi1D) {
				   if (lft1.dir != rgt2.dir) {
					   double q = (lft1.dir - rgt2.dir) * dk;
					   double qa = std::abs(q);
					   result = 2.0 * Nm1 * acutm1 * qa * gsl_sf_bessel_K1(qa * acut);
					   return result;
				   }
				   else {
					   return 2.0 * Nm1 * acutm2;
				   }
			   }
			   else {
				   return pm.intraMat[std::abs(lft1.dir - rgt2.dir)];
			   }
		   }
		   else {
			   if(!pm.quasi1D) {
				   if (lft1.dir != rgt2.dir) {
					   double q = (lft1.dir - rgt2.dir) * dk;
					   double qa = std::abs(q);
					   result =  pm.interMult * 2.0 * Nm1 * (Lzm1 * qa * gsl_sf_bessel_K1(qa * Lz) - qa * qa * gsl_sf_bessel_Kn(2, qa * Lz));
					   return result;
				   }
				   else {
					   return -pm.interMult * 2.0 * Nm1 * Lzm2;
				   }
			   }
			   else {
				   return pm.interMult * pm.interMat[std::abs(lft1.dir - rgt2.dir)];
			   }
		   }
	   }
	}
	
	return 0.0;
}
/*
void initIntralayerMat(double*& intraMat)
{
	intraMat = new double[NS];
	double q = 0.0;
	double value = 0.0;
	for (int qi = 0; qi < NS; ++qi) {
		value = 0.0;
		q = Lm1t * qi;
		value += (4 * dipConst / 3) * (1.0 / (dsz * dsy));
		if (qi == 0) {
			value -= 2.0 * dipConst / (dsy * (dsy + dsz));
		}
		else {
			value -= 2.0 * dipConst * std::exp(0.5 * q * q * dsz * dsz) * Inti(q);
		}
		intraMat[qi] = value;
	}
}

void initInterlayerMat(double*& interMat)
{
	interMat = new double[NS];
	double q = 0.0;
	double value = 0.0;
	for (int qi = 0; qi < NS; ++qi) {
		value = 0.0;
		q = Lm1t * qi;
		value += -(4.0 * dipConst / 3.0) * std::exp(-Lz * Lz / (2 * dsz * dsz)) / (2.0 * dsz * dsy);
		value += 2.0 * dipConst * std::exp(0.5 * q * q * dsy * dsy) * Intj(q);
		interMat[qi] = value;
		}
}


double IntI(double ky, void *paramsInt)
{
	intparams *intprm = (intparams *)paramsInt;
	double q = (intprm -> q);
	
	return sqrt(q * q + ky * ky) * std::exp(0.5 * ky * ky * (dsz * dsz - dsy * dsy)) * gsl_sf_erfc(sqrt(q * q + ky * ky) * dsz / sqrt(2.0));
}

double Inti(double q)
{
	gsl_integration_cquad_workspace* w = gsl_integration_cquad_workspace_alloc (10000);
	double min = 0.0;
	double max = 20.0;
	double eabs = 1.0e-7;
	double erel = 1.0e-7;
	double result = 0.0;
	double error;

	intparams intprm = {q};

	gsl_function F;
	
	F.function = &IntI;
	F.params = &intprm;

	gsl_integration_cquad(&F, min, max, eabs, erel, w, &result, &error, NULL);
	gsl_integration_cquad_workspace_free(w);
	return result;
}



double IntJ(double kz, void *paramsInt)
{
	intparams *intprm = (intparams *)paramsInt;
	double q = (intprm -> q);
	
	return kz * kz * cos(kz * Lz) * std::exp(0.5 * kz * kz * (dsy * dsy - dsz * dsz)) * gsl_sf_erfc(sqrt(q * q + kz * kz) * dsy / sqrt(2.0)) / sqrt(q * q + kz * kz);
}

double Intj(double q)
{
	gsl_integration_cquad_workspace* w = gsl_integration_cquad_workspace_alloc (10000);
	double min = 0.0;
	double max = 20.0;
	double eabs = 1.0e-7;
	double erel = 1.0e-7;
	double result = 0.0;
	double error;

	intparams intprm = {q};

	gsl_function F;
	
	F.function = &IntJ;
	F.params = &intprm;

	gsl_integration_cquad(&F, min, max, eabs, erel, w, &result, &error, NULL);
	gsl_integration_cquad_workspace_free(w);
	return result;
}
*/
void calcMat(SparseMat& baseHam, Count nimbsize, const State* qnummap, const Occup* qstates, const Param& pm)
{
	Count start =  baseHam.start;
	Count end = baseHam.end;
	Count valueSize  = 5 * nimbsize;
	MatType* baseHamMat = new MatType[valueSize];
	Count* baseHamia = new Count[end - start + 1];
	Count* baseHamja = new Count[valueSize];
	double matel;
	Count currbaseHamNum = 0;

	for (Count i = start; i < end; ++i) {
		for (Count j = i; j < nimbsize; ++j) {
			matel = MatElement(i, j, qnummap, qstates, pm);
			if (i == j) {
				baseHamia[i - start] = currbaseHamNum;
			}
			if (std::abs(matel) > thresh) {
				if (currbaseHamNum == valueSize) {
					valueSize *= 2;
					resize(baseHamMat, currbaseHamNum, valueSize);
					resize(baseHamja, currbaseHamNum, valueSize);
				}
				baseHamMat[currbaseHamNum] = matel;
				baseHamja[currbaseHamNum] = j;	
				++currbaseHamNum;
			}
		}
	}
	baseHamia[end - start] = currbaseHamNum;
	baseHam.mat = baseHamMat;
	baseHam.ia = baseHamia;
	baseHam.ja = baseHamja;
}

