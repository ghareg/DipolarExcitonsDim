#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include "basis.h"
#include "opermask.h"
#include "diag.h"
#include <thread>
#include <algorithm>
//#include "occupation.h"
#include "translation.h"
#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>


void sortEValues(eigen* eval, const MatType* evalues, const MatType* evecs, Count bSize, int neigsl);
void initStartEnd(SparseMat* baseHam, Count nimbsize);
void joinMat(SparseMat& MatFull, const SparseMat* MatArray, Count nimbsize, bool zeroBased);
void calcFullMat(SparseMat& baseHamFull, Count nimbsize, const State* qnummap, const Occup* qstates, const Param& pm);
template <typename arrayT>
void freeSparseMat(arrayT& Ham);

int main ()
{
	State* qnummap = nullptr;
	int qSize = 0;
	constrSingleParticleBasis(qnummap, qSize);

	SparseMat baseHamFull;
	GenMatProd op(baseHamFull);
/*	
	CMatType* redMat = new CMatType[(2 * NS - 1)];
	CMatType* redMatFourier = new CMatType[NS];
	CMatType* pairCorr = new CMatType[(2 * NS - 1)];
	CMatType* structFact = new CMatType[NS];
	int lyr1 = 0;
	int lyr2 = 1;
*/
	Count nimbsize = 0;
	Occup* qstates = nullptr;
	int NpL[NLayer] = {4, 0};
	for (int n = 1; n <= 10; ++n) {
		NpL[1] = n;
//		std::cout << NpL[0] << '\t' << NpL[1] << '\t' << lyr1 << '\t' << lyr2 << std::endl;
		
		std::cout << NpL[0] << '\t' << NpL[1] << std::endl;
		constrBasis(qstates, nimbsize, qnummap, qSize, NpL);
		//std::cout << nimbsize << std::endl;

		int neigsl = std::min(neigs, nimbsize / 4);
		MatType* evalues = new MatType[neigsl];
		MatType* evectors = new MatType[neigsl * nimbsize];
		MatType* evaluesWI = new MatType[neigsl];
		MatType* evectorsWI = new MatType[neigsl * nimbsize];

		MatType* momMat = new MatType[neigsl];

		eigen* eval = new eigen[neigsl];
		eigen* evalWI = new eigen[neigsl];

		Param pm;
		pm.tHop = tHop;
		pm.interMult = 1.0;
		pm.Ne = NpL[0] + n;
		pm.NpL[0] = NpL[0];
		pm.NpL[1] = NpL[1];
		pm.rhom1[0] =  1.0 * NS / NpL[0];
		pm.rhom1[1] =  1.0 * NS / NpL[1];
		pm.crossPair = false;
		pm.qSize = qSize;

//		initIntralayerMat(pm.intraMat);
//		initInterlayerMat(pm.interMat);
		pm.quasi1D = false;

		calcFullMat(baseHamFull, nimbsize, qnummap, qstates, pm);

		calcEValues(baseHamFull, op, evalues, evectors, neigsl);
		sortEValues(eval, evalues, evectors, nimbsize, neigsl);
		for (int i = 0; i < neigsl; ++i) {
			std::cout << '\t' << eval[i].value;
		}
		std::cout << std::endl;

		calcMomentumEigenvalues(momMat, neigsl, eval, qstates, nimbsize, qnummap, pm);
		for (int i = 0; i < neigsl; ++i) {
			std::cout << '\t' << momMat[i];
		}
		std::cout << std::endl;


/*		
		for (int indCM = 0; indCM < 0; NS; ++indCM) {
			pm.indCM = indCM;
		//	calcReducedDensMat(redMat, lyr1, lyr2, eval, qstates, nimbsize, qnummap, pm);
		//	printOccup(redMat, lyr1, lyr2);
		//	calcReducedDensMatFourier(redMatFourier, lyr1, lyr2, eval, qstates, nimbsize, qnummap, pm);
		//	printOccupFourier(redMatFourier, lyr1, lyr2);
			calcPairCorrelation(pairCorr, lyr1, lyr2, eval, qstates, nimbsize, qnummap, pm);
			printOccup(pairCorr, lyr1, lyr2);
		//	calcStructureFactor(structFact, lyr1, lyr2, eval, qstates, nimbsize, qnummap, pm);
		//	printOccupFourier(structFact, lyr1, lyr2);
			std::cout << std::endl;
		}*/
		freeSparseMat(baseHamFull);


		pm.interMult = 0.0;
		calcFullMat(baseHamFull, nimbsize, qnummap, qstates, pm);
		calcEValues(baseHamFull, op, evaluesWI, evectorsWI, neigsl);
		sortEValues(evalWI, evaluesWI, evectorsWI, nimbsize, neigsl);
		for (int i = 0; i < neigsl; ++i) {
			std::cout << '\t' << evalWI[i].value;
		}
		std::cout << std::endl;

		calcMomentumEigenvalues(momMat, neigsl, evalWI, qstates, nimbsize, qnummap, pm);
		for (int i = 0; i < neigsl; ++i) {
			std::cout << '\t' << momMat[i];
		}
		std::cout << std::endl;

/*		for (int indCM = 0; indCM < 0; NS; ++indCM) {
			pm.indCM = indCM;
		//	calcReducedDensMat(redMat, lyr1, lyr2, evalWI, qstates, nimbsize, qnummap, pm);
		//	printOccup(redMat, lyr1, lyr2);
		//	calcReducedDensMatFourier(redMatFourier, lyr1, lyr2, evalWI, qstates, nimbsize, qnummap, pm);
		//	printOccupFourier(redMatFourier, lyr1, lyr2);
			calcPairCorrelation(pairCorr, lyr1, lyr2, evalWI, qstates, nimbsize, qnummap, pm);
			printOccup(pairCorr, lyr1, lyr2);
		//	calcStructureFactor(structFact, lyr1, lyr2, evalWI, qstates, nimbsize, qnummap, pm);
		//	printOccupFourier(structFact, lyr1, lyr2);
			std::cout << std::endl;
		}
*/
		std::cout << eval[0].value - evalWI[0].value << std::endl;
/*
		delete[] pm.intraMat;
		delete[] pm.interMat;
*/

		delete[] qstates;
		qstates = nullptr;
		delete[] evalues;
		delete[] evectors;
		delete[] evaluesWI;
		delete[] evectorsWI;
		delete[] eval;
		delete[] evalWI;
		delete[] momMat;
		freeSparseMat(baseHamFull);
	}
/*
	delete[] redMat;
	delete[] redMatFourier;
	delete[] pairCorr;
	delete[] structFact;*/
	delete[] qnummap;
	return 0;
}

void sortEValues(eigen* eval, const MatType* evalues, const MatType* evecs, Count bSize, int neigsl)
{
	for (int i = 0; i < neigsl; ++i) {
		eval[i].value = evalues[i];
		eval[i].vector = &evecs[i * bSize];
	}
	std::sort(&eval[0], &eval[0] + neigsl);
}

void initStartEnd(SparseMat* baseHam, Count nimbsize)
{
	Count end0 = (Count) (nimbsize + 0.5 - 0.5 * sqrt((2.0 * nimbsize + 1.0) * (2.0 * nimbsize + 1.0) - 
			   4.0 * (1.0 * nimbsize * nimbsize + nimbsize) / nCore));
	baseHam[0].start = 0;
	baseHam[0].end = end0;
	for (int i = 1; i < nCore - 1; ++i) {
		Count oldstart = baseHam[i - 1].start;
		Count oldend = baseHam[i - 1].end;
		baseHam[i].start = oldend;
		Count endi = (Count) (nimbsize + 0.5 - 0.5 * sqrt((2.0 * nimbsize + 1.0) * (2.0 * nimbsize + 1.0) -
				16.0 * nimbsize * oldend - 8.0 * oldend + 8.0 * oldend * oldend + 8.0 * nimbsize * 
				oldstart - 4.0 * oldstart * oldstart + 4.0 * oldstart));
		baseHam[i].end = endi;
	}
	if (nCore > 1) {
		baseHam[nCore - 1].start = baseHam[nCore - 2].end;
	}
	baseHam[nCore - 1].end = nimbsize;
}


void joinMat(SparseMat& MatFull, const SparseMat* MatArray, Count bSize, bool zeroBased)
{
	Count valuesSize = 0;
	for (int i = 0; i < nCore; ++i) {
		valuesSize += MatArray[i].ia[MatArray[i].end - MatArray[i].start];
	}
	MatFull.mat = new MatType[valuesSize];
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

void calcFullMat(SparseMat& baseHamFull, Count nimbsize, const State* qnummap, const Occup* qstates, const Param& pm)
{
	SparseMat baseHam[nCore];
	initStartEnd(baseHam, nimbsize);
	std::thread t[nCore - 1];
	for (int i = 0; i < nCore - 1; ++i) {
		t[i] = std::thread(calcMat, std::ref(baseHam[i]), nimbsize, qnummap, qstates, std::ref(pm));
	}
			
	calcMat(baseHam[nCore - 1], nimbsize, qnummap, qstates, pm);

	for (int i = 0; i < nCore - 1; ++i) {
		t[i].join();
	}

	joinMat(baseHamFull, baseHam, nimbsize, false);
}

template <typename arrayT>
void freeSparseMat(arrayT& Ham)
{
	delete[] Ham.mat;
	delete[] Ham.ia;
	delete[] Ham.ja;
	Ham.mat = nullptr;
	Ham.ia = nullptr;
	Ham.ja = nullptr;
}
