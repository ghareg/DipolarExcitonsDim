#ifndef TRANSLATION_H_
#define TRANSLATION_H_

#include "basis.h"

void calcMomentumEigenvalues(MatType* momMat, int nsize, const eigen* eval, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm);

#endif
