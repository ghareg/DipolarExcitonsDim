#ifndef BASIS_H_
#define BASIS_H_
#include <stdint.h>
#include <complex>
#include "constants.h"

typedef int8_t Occup;
typedef int Count;
typedef int Index;
typedef double MatType;
typedef std::complex<MatType> CMatType;

const CMatType I(0.0, 1.0);

struct State
{
	State(): layer(0), dir(0) {}
	State(int lay, int pos) : layer(lay), dir(pos) {}
	Occup layer;
	Occup dir;
};

struct SparseMat
{
	SparseMat(): mat(nullptr), ia(nullptr), ja(nullptr), start(0), end(0) {}
	~SparseMat() {delete[] mat; delete[] ia; delete[] ja;}

	MatType* mat;
	Count* ia;
	Count* ja;
	int start;
	int end;
};

struct CSparseMat
{
	CSparseMat(): mat(nullptr), ia(nullptr), ja(nullptr), start(0), end(0) {}
	~CSparseMat() {delete[] mat; delete[] ia; delete[] ja;}

	CMatType* mat;
	Count* ia;
	Count* ja;
	int start;
	int end;
};

struct eigen {
	MatType value;
	const MatType *vector;
	
	bool operator<(eigen const &other) const {
		return value < other.value;
	}
};

struct Param
{
	double tHop;
	double interMult;
	int Ne;
	int NpL[NLayer];
	int qSize;
	double rhom1[NLayer];
	bool crossPair;
	int indCM;
	
	double* intraMat;
	double* interMat;
	bool quasi1D;

};

double factorial(int n);
void constrSingleParticleBasis(State*& qnummap, Count& qSize);
void constrBasis(Occup*& qstates, Count& stateNum, const State* qnummap, int qSize, const int* NpL);
void updateMBState(Occup* qstates, Count stateNum, Occup* occup, int Ne, int qSize);
void updateOccup(Occup* occup, Occup* layerOccup, const State* qnummap, int qSize, int Ne);
template <typename arrayT>
void resize(arrayT*& mat, Count currSize, Count newSize)
{
	arrayT* oldmat = mat;
	mat = new arrayT[newSize];
	for (Count i = 0; i < currSize; ++i) {
		mat[i] = oldmat[i];
	}
	delete[] oldmat;
}

#endif
