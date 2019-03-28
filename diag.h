#ifndef DIAG_H_
#define DIAG_H_
#include "basis.h"
#include <mkl_types.h>
#undef MKL_Complex16
#define MKL_Complex16 std::complex<double>
#include <mkl.h>

class GenMatProd
{
public:
	GenMatProd(const SparseMat& baseHam): baseHam_(baseHam) {}

	Count rows() {return baseHam_.end - baseHam_.start;}
	Count cols() {return baseHam_.end - baseHam_.start;}
	void perform_op(const MatType* xIn, MatType* yOut);
private:
	const SparseMat& baseHam_;
};

void calcEValues(const SparseMat& baseHam, GenMatProd& op, MatType* evalues, MatType* evecs, int neigs);

#endif
