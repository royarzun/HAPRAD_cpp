#ifndef TSEMIINCLUSIVEMODEL_H
#define TSEMIINCLUSIVEMODEL_H

#include "TROOT.h"

void SemiInclusiveModel(Double_t q2, Double_t X,
                        Double_t Y, Double_t Z,
                        Double_t pt2, Double_t mx2,
                        Double_t pl, Double_t& H1z,
                        Double_t& H2z, Double_t& H3z,
                        Double_t& H4z);

Double_t h3(Double_t X, Double_t q2, Double_t Z);

Double_t h4(Double_t X, Double_t q2, Double_t Z);

#endif


