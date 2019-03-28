#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const double Pi=3.1415926;
const int NLayer = 2;
const int NS = 8;
const double Lz = 0.5;
const double acut = 0.4;
const double acutm1 = 1.0 / acut;
const double acutm2 = 1.0 / (acut * acut);
const double Lzm1 = 1.0 / Lz;
const double Lzm2 = 1.0 / (Lz * Lz);
const double tHop = 0.0;

const double Pim1 = 1.0 / sqrt(Pi);
const double Nm1 = 1.0 / NS;

const double k0 = -Pi;
const double dk = 2.0 * Pi / NS;

const double dsz = 0.5 * Lz;
const double dsy = 0.5 * Lz;

const int neigs = 40;
const double thresh = 1E-6;
const int nCore = 4;

#endif
