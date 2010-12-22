#include "TKinematicalVariables.h"


TKinematicalVariables::TKinematicalVariables()
  : fX(0.), fY(0.), fZ(0.), fT(0.), fPhiH(0.), fE(0.)
{
    // Do nothing
}



TKinematicalVariables::~TKinematicalVariables()
{
    // Do nothing
}



TKinematicalVariables::TKinematicalVariables(Double_t x, Double_t Q2, Double_t z,
                                             Double_t p_t, Double_t phi, Double_t E)
  : fX(x), fY(Q2), fZ(z), fT(p_t), fPhiH(phi), fE(E)
{
    // Do nothing
}



void TKinematicalVariables::SetAll(Double_t x, Double_t y, Double_t z,
                                   Double_t t, Double_t phi, Double_t E)
{
    fX    = x;
    fY    = y;
    fZ    = z;
    fT    = t;
    fPhiH = phi;
    fE    = E;
}
