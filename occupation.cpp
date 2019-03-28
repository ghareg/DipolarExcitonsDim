#include "occupation.h"
#include <iostream>
#include "opermask.h"
#include <thread>


void calcReducedDensMat(CMatType* redMat, int lyr1, int lyr2, const eigen* eval, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm)
{
	int indc = 0;
	CSparseMat A;
	CMatType oc = 0.0;
	for (int ind = -NS + 1; ind < NS; ++ind) {
		ConstructFullInitMatrixRedDense(A, lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind, false);
		CalcProd(A, eval, nimbsize, oc, false);
		redMat[indc] = oc;
		++indc;
	}
}

void ConstructFullInitMatrixRedDense(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind, bool zeroBased)
{
	CSparseMat AHam[nCore];
	for (int i = 0; i < nCore - 1; ++i) {
		AHam[i].start = nimbsize * i / nCore;
		AHam[i].end = nimbsize * (i + 1) / nCore;
	}
	AHam[nCore - 1].start = nimbsize * (nCore - 1) / nCore;
	AHam[nCore - 1].end = nimbsize;

	std::thread t[nCore - 1];
	for (int i = 0; i < nCore - 1; ++i) {
		t[i] = std::thread(ConstructInitMatrixRedDense, std::ref(AHam[i]), lyr1, lyr2, qstates, nimbsize, qnummap, std::ref(pm), ind);
	}
			
	ConstructInitMatrixRedDense(AHam[nCore - 1], lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind);

	for (int i = 0; i < nCore - 1; ++i) {
		t[i].join();
	}
	
	joinMat(A, AHam, nimbsize, zeroBased);
}

void ConstructInitMatrixRedDense(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind)
{
	Count start =  A.start;
	Count end = A.end;
	Count valueSize  = 5 * nimbsize;
	A.mat = new CMatType[valueSize];
	A.ia = new Count[end - start + 1];
	A.ja = new Count[valueSize];
	Count currNum = 0;

	CMatType matel = CMatType(0.0, 0.0);
	for (Count row = start; row < end; ++row) {
		A.ia[row - start] = currNum;
		for (Count col = 0; col < nimbsize; ++col) {
			matel = RedDenseMatElement(row, col, lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind);
			if (std::abs(matel) > thresh) {
				if (currNum == valueSize) {
					valueSize *= 2;
					resize(A.mat, currNum, valueSize);
					resize(A.ja, currNum, valueSize);
				}
				A.mat[currNum] = matel;
				A.ja[currNum] = col;	
				++currNum;
			}
		}
	}
	A.ia[end - start] = currNum;
}

CMatType RedDenseMatElement(Count lftNum, Count rgtNum, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind)
{
	Index lft[2];
	Index rgt[2];
	int diffNum = constrDif(qstates, lftNum, rgtNum, lft, rgt, pm);
	if (diffNum == -1 || diffNum == 2) {
	   return 0.0;
	}

	double sum = 0.0;

	double sumCM = 0.0;
	CMatType matel(0.0, 0.0);
	if (lyr1 == lyr2 && diffNum == 0) {
		Occup const* lftOccup = &qstates[lftNum * pm.qSize];
		for (int i = 0; i < pm.qSize; ++i) {
			if (lftOccup[i] != 0 && qnummap[i].layer == lyr1) {
				sum = (qnummap[i].dir - 0.5 * NS) * ind * NSm1;
				matel += Nm1 * lftOccup[i] * std::exp(sum * I);
			}
		}
	}
	else if (diffNum == 1) {
		if (qnummap[lft[0]].layer == lyr1 && qnummap[rgt[0]].layer == lyr2) {
			if (qnummap[lft[0]].dir == qnummap[rgt[0]].dir) {
				sum = 0.5 * (qnummap[lft[0]].dir + qnummap[rgt[0]].dir - 0.5 * NS - 0.5 * NS) * ind * NSm1;
/*
			sumCM = (qnummap[lft[0]].dir - qnummap[rgt[0]].dir) * pm.indCM * NSm1;
			matel = Nm1 * sqrt(qstates[lftNum * pm.qSize + lft[0]] * qstates[rgtNum * pm.qSize + rgt[0]]) * std::exp(sum * I) * std::exp(sumCM * I);
*/			
				matel = Nm1 * sqrt(qstates[lftNum * pm.qSize + lft[0]] * qstates[rgtNum * pm.qSize + rgt[0]]) * std::exp(sum * I);
			}
		}
	}

	return matel;
}

void calcReducedDensMatFourier(CMatType* redMat, int lyr1, int lyr2, const eigen* eval, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm)
{
	int indc = 0;
	CSparseMat A;
	CMatType oc = 0.0;
	for (int ind = 0; ind < NS; ++ind) {
		ConstructFullInitMatrixRedDenseFourier(A, lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind, false);
		CalcProd(A, eval, nimbsize, oc, false);
		redMat[indc] = oc;
		++indc;
	}
}

void ConstructFullInitMatrixRedDenseFourier(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind, bool zeroBased)
{
	CSparseMat AHam[nCore];
	for (int i = 0; i < nCore - 1; ++i) {
		AHam[i].start = nimbsize * i / nCore;
		AHam[i].end = nimbsize * (i + 1) / nCore;
	}
	AHam[nCore - 1].start = nimbsize * (nCore - 1) / nCore;
	AHam[nCore - 1].end = nimbsize;

	std::thread t[nCore - 1];
	for (int i = 0; i < nCore - 1; ++i) {
		t[i] = std::thread(ConstructInitMatrixRedDenseFourier, std::ref(AHam[i]), lyr1, lyr2, qstates, nimbsize, qnummap, std::ref(pm), ind);
	}
			
	ConstructInitMatrixRedDenseFourier(AHam[nCore - 1], lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind);

	for (int i = 0; i < nCore - 1; ++i) {
		t[i].join();
	}
	
	joinMat(A, AHam, nimbsize, zeroBased);
}

void ConstructInitMatrixRedDenseFourier(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind)
{
	Count start =  A.start;
	Count end = A.end;
	Count valueSize  = 5 * nimbsize;
	A.mat = new CMatType[valueSize];
	A.ia = new Count[end - start + 1];
	A.ja = new Count[valueSize];
	Count currNum = 0;

	CMatType matel = CMatType(0.0, 0.0);
	for (Count row = start; row < end; ++row) {
		A.ia[row - start] = currNum;
		for (Count col = 0; col < nimbsize; ++col) {
			matel = RedDenseMatElementFourier(row, col, lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind);
			if (std::abs(matel) > thresh) {
				if (currNum == valueSize) {
					valueSize *= 2;
					resize(A.mat, currNum, valueSize);
					resize(A.ja, currNum, valueSize);
				}
				A.mat[currNum] = matel;
				A.ja[currNum] = col;
				++currNum;
			}
		}
	}
	A.ia[end - start] = currNum;
}

CMatType RedDenseMatElementFourier(Count lftNum, Count rgtNum, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind)
{
	Index lft[2];
	Index rgt[2];
	int diffNum = constrDif(qstates, lftNum, rgtNum, lft, rgt, pm);
	if (diffNum == -1 || diffNum == 2) {
	   return 0.0;
	}

	CMatType matel(0.0, 0.0);
	if (lyr1 == lyr2 && diffNum == 0) {
		Occup const* lftOccup = &qstates[lftNum * pm.qSize];
		for (int i = 0; i < pm.qSize; ++i) {
			if (lftOccup[i] != 0 && qnummap[i].layer == lyr1) {
				if (qnummap[i].dir == ind) {
					matel += lftOccup[i];
				}
			}
		}
	}
	else if (diffNum == 1) {
		if (qnummap[lft[0]].layer == lyr1 && qnummap[rgt[0]].layer == lyr2) {
			if (qnummap[lft[0]].dir == qnummap[rgt[0]].dir &&
					qnummap[lft[0]].dir == ind) {
				matel = sqrt(qstates[lftNum * pm.qSize + lft[0]] * qstates[rgtNum * pm.qSize + rgt[0]]);
			}
		}
	}

	return matel;
}

void calcPairCorrelation(CMatType* pairCorr, int lyr1, int lyr2, const eigen* eval, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm)
{
	int indc = 0;
	CSparseMat A;
	CMatType oc = 0.0;
	for (int ind = -NS + 1; ind < NS; ++ind) {
		ConstructFullPairCorr(A, lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind, false);
		CalcProd(A, eval, nimbsize, oc, false);
		pairCorr[indc] = oc;
		++indc;
	}
}

void ConstructFullPairCorr(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind, bool zeroBased)
{
	CSparseMat AHam[nCore];
	for (int i = 0; i < nCore - 1; ++i) {
		AHam[i].start = nimbsize * i / nCore;
		AHam[i].end = nimbsize * (i + 1) / nCore;
	}
	AHam[nCore - 1].start = nimbsize * (nCore - 1) / nCore;
	AHam[nCore - 1].end = nimbsize;

	std::thread t[nCore - 1];
	for (int i = 0; i < nCore - 1; ++i) {
		t[i] = std::thread(ConstructInitMatrixPairCorr, std::ref(AHam[i]), lyr1, lyr2, qstates, nimbsize, qnummap, std::ref(pm), ind);
	}
			
	ConstructInitMatrixPairCorr(AHam[nCore - 1], lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind);

	for (int i = 0; i < nCore - 1; ++i) {
		t[i].join();
	}
	
	joinMat(A, AHam, nimbsize, zeroBased);
}

void ConstructInitMatrixPairCorr(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind)
{
	Count start =  A.start;
	Count end = A.end;
	Count valueSize  = 5 * nimbsize;
	A.mat = new CMatType[valueSize];
	A.ia = new Count[end - start + 1];
	A.ja = new Count[valueSize];
	Count currNum = 0;

	CMatType matel = CMatType(0.0, 0.0);
	for (Count row = start; row < end; ++row) {
		A.ia[row - start] = currNum;
		for (Count col = 0; col < nimbsize; ++col) {
			matel = PairCorrMatElement(row, col, lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind);
			if (std::abs(matel) > thresh) {
				if (currNum == valueSize) {
					valueSize *= 2;
					resize(A.mat, currNum, valueSize);
					resize(A.ja, currNum, valueSize);
				}
				A.mat[currNum] = matel;
				A.ja[currNum] = col;	
				++currNum;
			}
		}
	}
	A.ia[end - start] = currNum;
}

CMatType PairCorrMatElement(Count lftNum, Count rgtNum, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind)
{
	Index lft[2];
	Index rgt[2];
	int diffNum = constrDif(qstates, lftNum, rgtNum, lft, rgt, pm);
	if (diffNum == -1) {
	   return 0.0;
	}

	CMatType matel(0.0, 0.0);
	Occup const* lftOccup = &qstates[lftNum * pm.qSize];
	Occup const* rgtOccup = &qstates[rgtNum * pm.qSize];
	if (diffNum == 0) {
		for (int i = 0; i < pm.qSize; ++i) {
			if (lftOccup[i] > 0) {
				for (int j = i + 1; j < pm.qSize; ++j) {
					if (lftOccup[j] > 0) {
						matel += lftOccup[i] * lftOccup[j] * 1.0 * (PairTwoMatEl(qnummap[i], qnummap[j], qnummap[i], qnummap[j], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[i], qnummap[j], qnummap[j], qnummap[i], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[j], qnummap[i], qnummap[i], qnummap[j], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[j], qnummap[i], qnummap[j], qnummap[i], lyr1, lyr2, ind, pm));
					}
				}
				if (lftOccup[i] > 1) {
					matel += lftOccup[i] * (lftOccup[i] - 1) * 1.0 * PairTwoMatEl(qnummap[i], qnummap[i], qnummap[i], qnummap[i], lyr1, lyr2, ind, pm);
				}
			}
		}
	}
	else if (diffNum == 1) {
		for (int i = 0; i < pm.qSize; ++i) {
			if (i == lft[0]) {
				if (lftOccup[i] > 1 && rgtOccup[i] > 0) {
					matel += sqrt(lftOccup[i] * (lftOccup[i] - 1) * rgtOccup[i] * rgtOccup[rgt[0]]) * (PairTwoMatEl(qnummap[i], qnummap[i], qnummap[i], qnummap[rgt[0]], lyr1, lyr2, ind, pm) +  PairTwoMatEl(qnummap[i], qnummap[i], qnummap[rgt[0]], qnummap[i], lyr1, lyr2, ind, pm));
				}
			}
			else if (i == rgt[0]) {
				if (lftOccup[i] > 0 && rgtOccup[i] > 1) {
					matel += sqrt(lftOccup[i] * lftOccup[lft[0]] * rgtOccup[i] * (rgtOccup[i] - 1)) * (PairTwoMatEl(qnummap[lft[0]], qnummap[i], qnummap[i], qnummap[i], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[i], qnummap[lft[0]], qnummap[i], qnummap[i], lyr1, lyr2, ind, pm));
				}
			}
			else {
				if (lftOccup[i] > 0 && rgtOccup[i] > 0) {
					matel += sqrt(lftOccup[i] * lftOccup[lft[0]] * rgtOccup[i] * rgtOccup[rgt[0]]) * (PairTwoMatEl(qnummap[lft[0]], qnummap[i], qnummap[i], qnummap[rgt[0]], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[lft[0]], qnummap[i], qnummap[rgt[0]], qnummap[i], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[i], qnummap[lft[0]], qnummap[i], qnummap[rgt[0]], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[i], qnummap[lft[0]], qnummap[rgt[0]], qnummap[i], lyr1, lyr2, ind, pm));
				}
			}
		}
	}
	else {
		if (lft[0] == lft[1] && rgt[0] == rgt[1]) {
			matel += sqrt(lftOccup[lft[0]] * (lftOccup[lft[0]] - 1) * rgtOccup[rgt[0]] * (rgtOccup[rgt[0]] - 1)) * PairTwoMatEl(qnummap[lft[0]], qnummap[lft[0]], qnummap[rgt[0]], qnummap[rgt[0]], lyr1, lyr2, ind, pm);
		}
		else if (lft[0] == lft[1]) {
			matel += sqrt(lftOccup[lft[0]] * (lftOccup[lft[0]] - 1) * rgtOccup[rgt[0]] * rgtOccup[rgt[1]]) * (PairTwoMatEl(qnummap[lft[0]], qnummap[lft[0]], qnummap[rgt[0]], qnummap[rgt[1]], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[lft[0]], qnummap[lft[0]], qnummap[rgt[1]], qnummap[rgt[0]], lyr1, lyr2, ind, pm)); 
		}
		else if (rgt[0] == rgt[1]) {
			matel += sqrt(lftOccup[lft[0]] * lftOccup[lft[1]] * rgtOccup[rgt[0]] * (rgtOccup[rgt[0]] - 1)) * (PairTwoMatEl(qnummap[lft[0]], qnummap[lft[1]], qnummap[rgt[0]], qnummap[rgt[0]], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[lft[1]], qnummap[lft[0]], qnummap[rgt[0]], qnummap[rgt[0]], lyr1, lyr2, ind, pm));  
		}
		else {
			matel += sqrt(lftOccup[lft[0]] * lftOccup[lft[1]] * rgtOccup[rgt[0]] * rgtOccup[rgt[1]]) * (PairTwoMatEl(qnummap[lft[0]], qnummap[lft[1]], qnummap[rgt[0]], qnummap[rgt[1]], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[lft[0]], qnummap[lft[1]], qnummap[rgt[1]], qnummap[rgt[0]], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[lft[1]], qnummap[lft[0]], qnummap[rgt[0]], qnummap[rgt[1]], lyr1, lyr2, ind, pm) + PairTwoMatEl(qnummap[lft[1]], qnummap[lft[0]], qnummap[rgt[1]], qnummap[rgt[0]], lyr1, lyr2, ind, pm));  
		}
	}

	return matel;
}
	
CMatType PairTwoMatEl(const State& lft1, const State& lft2, const State& rgt1, const State& rgt2, int lyr1, int lyr2, int ind, const Param& pm)
{
	double sum = 0.0;

	//double sumCM = 0.0;
	if (!pm.crossPair) {
		if (lft1.layer == rgt2.layer && lft2.layer == rgt1.layer &&
			   lft1.layer == lyr1 && lft2.layer == lyr2) {
			if (lft1.dir + lft2.dir == rgt1.dir + rgt2.dir) {
				sum = 0.5 * (lft1.dir - lft2.dir + rgt1.dir - rgt2.dir) * ind * NSm1;
			
			/*sumCM = (lft1.dir + lft2.dir - rgt1.dir - rgt2.dir) * pm.indCM * NSm1;

			return Nm1 * Nm1 * pm.rhom1[lyr1] * pm.rhom1[lyr2] * std::exp(sum * I) * std::exp(sumCM * I);
			*/
				return Nm1 * Nm1 * pm.rhom1[lyr1] * pm.rhom1[lyr2] * std::exp(sum * I);
			}
		}
	}
	else {
		if (lft1.layer == rgt1.layer && lft2.layer == rgt2.layer &&
			   lft1.layer == lyr1 && lft2.layer == lyr2) {
			if (lft1.dir + lft2.dir == rgt1.dir + rgt2.dir) {
				sum = 0.5 * (lft1.dir - lft2.dir + rgt1.dir - rgt2.dir) * ind * NSm1;
			/*
			sumCM = (lft1.dir + lft2.dir - rgt1.dir - rgt2.dir) * pm.indCM * NSm1;

			return Nm1 * Nm1 * pm.rhom1[lyr1] * pm.rhom1[lyr2] * std::exp(sum * I) * std::exp(sumCM * I);
			*/

				return Nm1 * Nm1 * pm.rhom1[lyr1] * pm.rhom1[lyr2] * std::exp(sum * I);
			}
		}
	}
	return 0.0;
}

void calcStructureFactor(CMatType* structFact, int lyr1, int lyr2, const eigen* eval, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm)
{
	int indc = 0;
	CSparseMat A;
	CMatType oc = 0.0;
	double sum = 0.0;
	for (int ind = 0; ind < NS; ++ind) {
		ConstructFullStructFact(A, lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind, false);
		CalcProd(A, eval, nimbsize, oc, false);
		sum = 0.0;
		if (lyr1 == lyr2 || pm.crossPair) {
			sum += 1.0 * pm.NpL[lyr1] / pm.Ne;
		}
		if (ind == 0) {
			sum += -1.0 * pm.NpL[lyr1] * pm.NpL[lyr2] / pm.Ne;
		}
		structFact[indc] = oc + sum;
		++indc;
	}
}

void ConstructFullStructFact(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind, bool zeroBased)
{
	CSparseMat AHam[nCore];
	for (int i = 0; i < nCore - 1; ++i) {
		AHam[i].start = nimbsize * i / nCore;
		AHam[i].end = nimbsize * (i + 1) / nCore;
	}
	AHam[nCore - 1].start = nimbsize * (nCore - 1) / nCore;
	AHam[nCore - 1].end = nimbsize;

	std::thread t[nCore - 1];
	for (int i = 0; i < nCore - 1; ++i) {
		t[i] = std::thread(ConstructInitMatrixStructFact, std::ref(AHam[i]), lyr1, lyr2, qstates, nimbsize, qnummap, std::ref(pm), ind);
	}
			
	ConstructInitMatrixStructFact(AHam[nCore - 1], lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind);

	for (int i = 0; i < nCore - 1; ++i) {
		t[i].join();
	}
	
	joinMat(A, AHam, nimbsize, zeroBased);
}

void ConstructInitMatrixStructFact(CSparseMat& A, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind)
{
	Count start =  A.start;
	Count end = A.end;
	Count valueSize  = 5 * nimbsize;
	A.mat = new CMatType[valueSize];
	A.ia = new Count[end - start + 1];
	A.ja = new Count[valueSize];
	Count currNum = 0;

	CMatType matel = CMatType(0.0, 0.0);
	for (Count row = start; row < end; ++row) {
		A.ia[row - start] = currNum;
		for (Count col = 0; col < nimbsize; ++col) {
			matel = StructFactMatElement(row, col, lyr1, lyr2, qstates, nimbsize, qnummap, pm, ind);
			if (std::abs(matel) > thresh) {
				if (currNum == valueSize) {
					valueSize *= 2;
					resize(A.mat, currNum, valueSize);
					resize(A.ja, currNum, valueSize);
				}
				A.mat[currNum] = matel;
				A.ja[currNum] = col;	
				++currNum;
			}
		}
	}
	A.ia[end - start] = currNum;
}

CMatType StructFactMatElement(Count lftNum, Count rgtNum, int lyr1, int lyr2, const Occup* qstates, Count nimbsize, const State* qnummap, const Param& pm, int ind)
{
	Index lft[2];
	Index rgt[2];
	int diffNum = constrDif(qstates, lftNum, rgtNum, lft, rgt, pm);
	if (diffNum == -1) {
	   return 0.0;
	}

	CMatType matel(0.0, 0.0);
	Occup const* lftOccup = &qstates[lftNum * pm.qSize];
	Occup const* rgtOccup = &qstates[rgtNum * pm.qSize];
	if (diffNum == 0) {
		for (int i = 0; i < pm.qSize; ++i) {
			if (lftOccup[i] > 0) {
				for (int j = i + 1; j < pm.qSize; ++j) {
					if (lftOccup[j] > 0) {
						matel += lftOccup[i] * lftOccup[j] * 1.0 * (StructTwoMatEl(qnummap[i], qnummap[j], qnummap[i], qnummap[j], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[i], qnummap[j], qnummap[j], qnummap[i], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[j], qnummap[i], qnummap[i], qnummap[j], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[j], qnummap[i], qnummap[j], qnummap[i], lyr1, lyr2, ind, pm));
					}
				}
				if (lftOccup[i] > 1) {
					matel += lftOccup[i] * (lftOccup[i] - 1) * 1.0 * StructTwoMatEl(qnummap[i], qnummap[i], qnummap[i], qnummap[i], lyr1, lyr2, ind, pm);
				}
			}
		}
	}
	else if (diffNum == 1) {
		for (int i = 0; i < pm.qSize; ++i) {
			if (i == lft[0]) {
				if (lftOccup[i] > 1 && rgtOccup[i] > 0) {
					matel += sqrt(lftOccup[i] * (lftOccup[i] - 1) * rgtOccup[i] * rgtOccup[rgt[0]]) * (StructTwoMatEl(qnummap[i], qnummap[i], qnummap[i], qnummap[rgt[0]], lyr1, lyr2, ind, pm) +  StructTwoMatEl(qnummap[i], qnummap[i], qnummap[rgt[0]], qnummap[i], lyr1, lyr2, ind, pm));
				}
			}
			else if (i == rgt[0]) {
				if (lftOccup[i] > 0 && rgtOccup[i] > 1) {
					matel += sqrt(lftOccup[i] * lftOccup[lft[0]] * rgtOccup[i] * (rgtOccup[i] - 1)) * (StructTwoMatEl(qnummap[lft[0]], qnummap[i], qnummap[i], qnummap[i], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[i], qnummap[lft[0]], qnummap[i], qnummap[i], lyr1, lyr2, ind, pm));
				}
			}
			else {
				if (lftOccup[i] > 0 && rgtOccup[i] > 0) {
					matel += sqrt(lftOccup[i] * lftOccup[lft[0]] * rgtOccup[i] * rgtOccup[rgt[0]]) * (StructTwoMatEl(qnummap[lft[0]], qnummap[i], qnummap[i], qnummap[rgt[0]], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[lft[0]], qnummap[i], qnummap[rgt[0]], qnummap[i], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[i], qnummap[lft[0]], qnummap[i], qnummap[rgt[0]], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[i], qnummap[lft[0]], qnummap[rgt[0]], qnummap[i], lyr1, lyr2, ind, pm));
				}
			}
		}
	}
	else {
		if (lft[0] == lft[1] && rgt[0] == rgt[1]) {
			matel += sqrt(lftOccup[lft[0]] * (lftOccup[lft[0]] - 1) * rgtOccup[rgt[0]] * (rgtOccup[rgt[0]] - 1)) * StructTwoMatEl(qnummap[lft[0]], qnummap[lft[0]], qnummap[rgt[0]], qnummap[rgt[0]], lyr1, lyr2, ind, pm);
		}
		else if (lft[0] == lft[1]) {
			matel += sqrt(lftOccup[lft[0]] * (lftOccup[lft[0]] - 1) * rgtOccup[rgt[0]] * rgtOccup[rgt[1]]) * (StructTwoMatEl(qnummap[lft[0]], qnummap[lft[0]], qnummap[rgt[0]], qnummap[rgt[1]], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[lft[0]], qnummap[lft[0]], qnummap[rgt[1]], qnummap[rgt[0]], lyr1, lyr2, ind, pm)); 
		}
		else if (rgt[0] == rgt[1]) {
			matel += sqrt(lftOccup[lft[0]] * lftOccup[lft[1]] * rgtOccup[rgt[0]] * (rgtOccup[rgt[0]] - 1)) * (StructTwoMatEl(qnummap[lft[0]], qnummap[lft[1]], qnummap[rgt[0]], qnummap[rgt[0]], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[lft[1]], qnummap[lft[0]], qnummap[rgt[0]], qnummap[rgt[0]], lyr1, lyr2, ind, pm));  
		}
		else {
			matel += sqrt(lftOccup[lft[0]] * lftOccup[lft[1]] * rgtOccup[rgt[0]] * rgtOccup[rgt[1]]) * (StructTwoMatEl(qnummap[lft[0]], qnummap[lft[1]], qnummap[rgt[0]], qnummap[rgt[1]], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[lft[0]], qnummap[lft[1]], qnummap[rgt[1]], qnummap[rgt[0]], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[lft[1]], qnummap[lft[0]], qnummap[rgt[0]], qnummap[rgt[1]], lyr1, lyr2, ind, pm) + StructTwoMatEl(qnummap[lft[1]], qnummap[lft[0]], qnummap[rgt[1]], qnummap[rgt[0]], lyr1, lyr2, ind, pm));  
		}
	}

	return matel;
}
	
CMatType StructTwoMatEl(const State& lft1, const State& lft2, const State& rgt1, const State& rgt2, int lyr1, int lyr2, int ind, const Param& pm)
{
	if (!pm.crossPair) {
		if (lft1.layer == rgt2.layer && lft2.layer == rgt1.layer &&
			   lft1.layer == lyr1 && lft2.layer == lyr2) {
			if (lft1.dir + lft2.dir == rgt1.dir + rgt2.dir) {
				if (lft1.dir - rgt2.dir == ind && rgt1.dir - lft2.dir == ind) {
					return 1.0 / pm.Ne;
				}
			}
		}
	}
	else {
		if (lft1.layer == rgt1.layer && lft2.layer == rgt2.layer &&
			   lft1.layer == lyr1 && lft2.layer == lyr2) {
			if (lft1.dir + lft2.dir == rgt1.dir + rgt2.dir) {
				if (lft1.dir - rgt2.dir == ind && rgt1.dir - lft2.dir == ind) {
					return 1 / pm.Ne;
				}
			}
		}
	}
	return 0.0;
}


void CalcProd(const CSparseMat& A, const eigen* eval, Count nimbsize, CMatType& oc, bool zeroBased)
{
	CMatType* evecs = new CMatType[nimbsize];
	CMatType* res = new CMatType[nimbsize];
	for (Count ind = 0; ind < nimbsize; ++ind) {
		evecs[ind] = eval[0].vector[ind];
	}
	char transa = 'N';
	CMatType ocr(0.0, 0.0);
	if (zeroBased) {
		mkl_cspblas_zcsrgemv (&transa, &nimbsize, A.mat, A.ia, A.ja, evecs, res);
	}
	else {
		mkl_zcsrgemv (&transa, &nimbsize, A.mat, A.ia, A.ja, evecs, res);
	}
	cblas_zdotc_sub (nimbsize, (void *)evecs, 1, (void *) res, 1, (void *) &ocr);
	oc = ocr;
	delete[] evecs;
	delete[] res;
}

void joinMat(CSparseMat& MatFull, const CSparseMat* MatArray, Count bSize, bool zeroBased)
{
	Count valuesSize = 0;
	for (int i = 0; i < nCore; ++i) {
		valuesSize += MatArray[i].ia[MatArray[i].end - MatArray[i].start];
	}
	MatFull.mat = new CMatType[valuesSize];
	MatFull.ja = new Count[valuesSize];
	MatFull.ia = new Count[bSize + 1];
	Count jadd = zeroBased ? 0 : 1;
	Count indCount = 0;
	Count iaIndCount = 0;
	Count currTot = zeroBased ? 0 : 1;
	for (int i = 0; i < nCore; ++i) {
		for (Count j = 0; j < MatArray[i].ia[MatArray[i].end - MatArray[i].start]; ++j) {
			MatFull.mat[indCount] = MatArray[i].mat[j];
			MatFull.ja[indCount] = MatArray[i].ja[j] + jadd;
			indCount++;
		}
		for (Count j = 0; j < MatArray[i].end - MatArray[i].start; ++j) {
			MatFull.ia[iaIndCount] = currTot + MatArray[i].ia[j];
			++iaIndCount;
		}
		currTot += MatArray[i].ia[MatArray[i].end - MatArray[i].start];
	}
	
	MatFull.ia[bSize] = currTot; 
	MatFull.start = 0;
	MatFull.end = bSize;
}

void diagtoFullMat(CSparseMat& MatFull, CSparseMat& MatUFull, bool zeroBased)
{
	const int add = zeroBased ? 0 : 1;
	const Count nimbsize = MatUFull.end - MatUFull.start;
	const Count USize = MatUFull.ia[nimbsize] - add;
	CSparseMat MatDFull;
	MatDFull.mat = new CMatType[USize];
	MatDFull.ja = new Count[USize];
	MatDFull.ia = new Count[nimbsize + 1];
	MKL_INT job[6];
	job[0] = 1;
	job[1] = add;
	job[2] = add;
	job[5] = 1;
	const MKL_INT n = nimbsize;
	MKL_INT info;
	mkl_zcsrcsc (job, &n , MatDFull.mat, MatDFull.ja , MatDFull.ia , MatUFull.mat, MatUFull.ja , MatUFull.ia , &info);
	for (Count i = 0; i < USize; ++i) {
		MatDFull.mat[i] = std::conj(MatDFull.mat[i]);
	}
	for (Count i = 0; i < nimbsize; ++i) {
		MatUFull.mat[MatUFull.ia[i] - add] = CMatType(0.0, 0.0);
	}
	Count Size = 2 * (USize - nimbsize) + nimbsize;
	MatFull.mat = new CMatType[Size];
	MatFull.ja = new Count[Size];
	MatFull.ia = new Count[nimbsize + 1];
	MatFull.start = MatUFull.start;
	MatFull.end = MatUFull.end;
	const char trans = 'N';
	MKL_INT request = 0;
	CMatType beta = CMatType(1.0, 0.0);
	const MKL_INT nzmax = Size;
	MKL_INT sort = 0;
	mkl_zcsradd (&trans, &request, &sort, &n , &n , MatUFull.mat, MatUFull.ja, MatUFull.ia, &beta, MatDFull.mat, MatDFull.ja, MatDFull.ia, MatFull.mat, MatFull.ja, MatFull.ia, &nzmax, &info);
}

void printOccup(const CMatType* occup, int lyr1, int lyr2)
{
	int indc = 0;
//	std::cout << lyr1 << "\t" << lyr2 << std::endl;
	for (int ind = -NS + 1; ind < NS; ++ind) {
		std::cout << ind << "\t" << occup[indc].real() << "\t" << occup[indc].imag() << std::endl;
		++indc;
	}
}

void printOccupFourier(const CMatType* occup, int lyr1, int lyr2)
{
	int indc = 0;
//	std::cout << lyr1 << "\t" << lyr2 << std::endl;
	for (int ind = 0; ind < NS; ++ind) {
		std::cout << ind << "\t" << occup[indc].real() << "\t" << occup[indc].imag() << std::endl;
		++indc;
	}
}
