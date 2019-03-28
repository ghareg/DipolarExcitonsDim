#ifndef OCCUPATION_H_
#define OCCUPATION_H_
#include "basis.h"
#include <complex>
#include <mkl_types.h>
#undef MKL_Complex16
#define MKL_Complex16 std::complex<double>
#include <mkl.h>

void calcReducedDensMat(CMatType* redMat, int lyr1, int lyr2, const eigen* eval, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm);
void ConstructFullInitMatrixRedDense(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind, bool zeroBased);
void ConstructInitMatrixRedDense(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind);
CMatType RedDenseMatElement(Count lftNum, Count rgtNum, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind);

void calcReducedDensMatFourier(CMatType* redMat, int lyr1, int lyr2, const eigen* eval, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm);
void ConstructFullInitMatrixRedDenseFourier(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind, bool zeroBased);
void ConstructInitMatrixRedDenseFourier(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind);
CMatType RedDenseMatElementFourier(Count lftNum, Count rgtNum, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind);


void calcPairCorrelation(CMatType* pairCorr, int lyr1, int lyr2, const eigen* eval, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm);
void ConstructFullPairCorr(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind, bool zeroBased);
void ConstructInitMatrixPairCorr(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind);
CMatType PairCorrMatElement(Count lftNum, Count rgtNum, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind);
CMatType PairTwoMatEl(const State& lft1, const State& lft2, const State& rgt1, const State& rgt2, int lyr1, int lyr2, int ind, const Param& pm);

void calcStructureFactor(CMatType* structFact, int lyr1, int lyr2, const eigen* eval, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm);
void ConstructFullStructFact(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind, bool zeroBased);
void ConstructInitMatrixStructFact(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind);
CMatType StructFactMatElement(Count lftNum, Count rgtNum, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind);
CMatType StructTwoMatEl(const State& lft1, const State& lft2, const State& rgt1, const State& rgt2, int lyr1, int lyr2, int ind, const Param& pm);

void CalcProd(const CSparseMat& A, const eigen* eval, Count nimbsize, CMatType& oc, bool zeroBased);
void joinMat(CSparseMat& MatFull, const CSparseMat* MatArray, Count bSize, bool zeroBased);
void diagtoFullMat(CSparseMat& MatFull, CSparseMat& MatUFull, bool zeroBased);
void printOccup(const CMatType* occup, int lyr1, int lyr2);
void printOccupFourier(const CMatType* occup, int lyr1, int lyr2);

#endif
