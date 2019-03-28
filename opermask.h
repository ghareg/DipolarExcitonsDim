#ifndef OPERMASK_H_
#define OPERMASK_H_

#include "basis.h"

int constrDif(const Occup* qstates, Count lftNum, Count rgtNum, Index* lft, Index* rgt, const Param& pm);
double MatElement(Count lftNum, Count rgtNum, const State* qnummap, const Occup* qstates, const Param& pm);
double OneBodyElement(const State& st, const Param& pm);
double TwoBodyElement(const State& lft1, const State& lft2, const State& rgt1, const State& rgt2, const Param& pm);
/*
void initIntralayerMat(double*& intraMat);
void initInterlayerMat(double*& interMat);
double IntI(double ky, void *paramsInt);
double Inti(double q);
double IntJ(double kz, void *paramsInt);
double Intj(double q);
struct intparams {double q;};
*/
void calcMat(SparseMat& baseHam, Count nimbsize, const State* qnummap, const Occup* qstates, const Param& pm);

#endif
