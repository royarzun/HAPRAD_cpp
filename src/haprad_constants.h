#ifndef HAPRAD_CONSTANTS_H
#define HAPRAD_CONSTANTS_H

#include "TMath.h"


static const Double_t pi = 3.1415926;
//static const Double_t amhh = 0.1395675;
static const Double_t amp = 0.938272;
static const Double_t ampi = 0.1395675;
//static const Double_t amlep[2] = {0.511000*TMath::Power(10,-3),0.10565};

static const Int_t i1[8] = {3,3,3,3,3,3,3,3};
static const Int_t i2[8] = {1,1,1,1,1,1,1,1};
static const Int_t isf20 = 4;

//from haprad: static constants8.inc
static const Double_t mp = 0.938272;
//static const Double_t me = 0.5110*TMath::Power(10,-3);
static const Double_t mpi = 0.1395675;
//static const Double_t alpha = 7.297*TMath::Power(10,-3);
static const Double_t mn = 0.93956536;
#endif
