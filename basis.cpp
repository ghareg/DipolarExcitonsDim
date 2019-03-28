#include "basis.h"

double factorial(int n)
{
	double sum = 1.0;
	for(int i = 1; i <=n; ++i) {
		sum *= i;
	}
	return sum;
}

void constrSingleParticleBasis(State*& qnummap, Count& qSize)
{
	int qMaxSize = NLayer * NS;
	qnummap = new State[qMaxSize];

	qSize = 0;
	for (int layer = 0; layer < NLayer; ++layer) {
		for (int pos = 0; pos < NS; ++pos) {
			qnummap[qSize++] = State(layer, pos);
		}
	}
}

void constrBasis(Occup*& qstates, Count& stateNum, const State* qnummap, int qSize, const int* NpL)
{
	int Ne = 0;
	Count qstMaxSize = 1;
	for (int layer = 0; layer < NLayer; ++layer) {
		qstMaxSize *= Count(factorial(NpL[layer] + NS - 1) / (factorial(NpL[layer]) * factorial(NS - 1)));
		Ne += NpL[layer];
	}
	qstates = new Occup[qstMaxSize * qSize];
	
	stateNum = 0;
	Occup layerOccup[NLayer] = {0};
	Occup occup[Ne];
	for (int i = 0; i < Ne; ++i) {
		occup[i] = 0;
	}
	layerOccup[qnummap[0].layer] = Ne;
	bool accState = true;
	while (occup[0] < qSize) {
		accState = true;
		for (int layer = 0; layer < NLayer; ++layer) {
			if (layerOccup[layer] != NpL[layer]) {
				accState = false;
				break;
			}
		}
		if (accState) {
			if (qstMaxSize == stateNum) {
				qstMaxSize = Count(1.2 * qstMaxSize);
				resize(qstates, stateNum * qSize, qstMaxSize * qSize);
			}
			updateMBState(qstates, stateNum, occup, Ne, qSize);
			++stateNum;
		}
		updateOccup(occup, layerOccup, qnummap, qSize, Ne);
	}
	resize(qstates, stateNum * qSize, stateNum * qSize);
}

void updateMBState(Occup* qstates, Count stateNum, Occup* occup, int Ne, int qSize)
{
	for (int i = 0; i < qSize; ++i) {
		qstates[stateNum * qSize + i] = 0;
	}
	int i = 0;
	int state = 0;
	int count = 0;
	while (i < Ne) {
		state = occup[i];
		count = 1;
		++i;
		while (i < Ne && state == occup[i]) {
			++i;
			++count;
		}
		qstates[stateNum * qSize + state] += count;
	}
}

void updateOccup(Occup* occup, Occup* layerOccup, const State* qnummap, int qSize, int Ne)
{
	int ind = Ne - 1;
	while (ind > 0 && occup[ind] == qSize - 1) {
		--ind;
	}
	if (ind < Ne - 1) {
		layerOccup[qnummap[qSize - 1].layer] -= Ne - 1 - ind;
	}
	
	--layerOccup[qnummap[occup[ind]].layer];
	++occup[ind];
	int initNum = occup[ind];
	
	for (int indCur = 1; indCur < Ne - ind; ++indCur) {
		occup[ind + indCur] = initNum;
	}
	if (initNum < qSize) {
		layerOccup[qnummap[initNum].layer] += Ne - ind;
	}
}
