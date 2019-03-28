#include "translation.h"

void calcMomentumEigenvalues(MatType* momMat, int nsize, const eigen* eval, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm)
{
//	double a = 1.0;
	CMatType* MultVector = new CMatType[nimbsize];
	CMatType mult = 1.0;
	for (Count ind = 0; ind < nimbsize; ++ind) {
		mult = 1.0;
		const Occup* curOccup = &qstates[ind * pm.qSize];
		for (int i = 0; i < pm.qSize; ++i) {
			if (curOccup[i] != 0) {
				mult *= std::exp(-curOccup[i] * (k0 + qnummap[i].dir * dk) * I);
			}
		}
		MultVector[ind] = mult;
	}
	CMatType result = 0.0;
	for (int i = 0; i < nsize; ++i) {
		result = 0.0;
		for (Count ind = 0; ind < nimbsize; ++ind) {
			result += MultVector[ind] * std::norm(eval[i].vector[ind]);
		}
		momMat[i] = atan2(-result.imag(), result.real());
	}

	delete[] MultVector;
}
