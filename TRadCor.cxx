#include "TRadCor.h"
#include "haprad_constants.h"
#include <cmath>


TRadCor::TRadCor()
: fEbeam(0),
  fX(0), fY(0), fZ(0), fPt(0), fPhi(0),
  fMaxMx2(0), fMx2(0),
  fSib(0), fSig(0), fDelta(0), fTail(0)
{
    // Default constructor
}


TRadCor::TRadCor(Double_t Ebeam, Double_t x, Double_t q2, Double_t z,
                 Double_t pt, Double_t phi, Double_t maxMx2)
: fEbeam(Ebeam),
  fX(x), fY(-q2), fZ(z), fPt(pt), fPhi(phi/kRadianDeg),
  fMaxMx2(maxMx2), fMx2(0),
  fSib(0), fSig(0), fDelta(0), fTail(0)
{
    // Normal constructor for a radiative correction object
    //
    // The Ebeam parameter is the energy of the beam, the x, q2, z, pt and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    Setup();
}



TRadCor::~TRadCor()
{
    // Default destructor
}



void TRadCor::SetParameters(Double_t Ebeam, Double_t x, Double_t q2, Double_t z,
                            Double_t pt, Double_t phi, Double_t maxMx2)
{
    // Set the values for the variables used to calculate the radiative
    // correction.
    //
    // The Ebeam parameter is the energy of the beam, the x, q2, z, pt and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    fEbeam = Ebeam;
    fX = x;
    fY = -q2;
    fZ = z;
    fPt = pt;
    fPhi = phi / kRadianDeg;
    fMaxMx2 = maxMx2;

    Setup();
}



void TRadCor::SetEbeam(Double_t Ebeam)
{
    fEbeam = Ebeam;
}



void TRadCor::SetX(Double_t x)
{
    fX = x;
    Setup();
}



void TRadCor::SetQ2(Double_t q2)
{
    fY = -q2;
    Setup();
}



void TRadCor::SetZ(Double_t z)
{
    fZ = z;
    Setup();
}



void TRadCor::SetPt(Double_t pt)
{
    fPt = pt;
    Setup();
}



void TRadCor::SetPhi(Double_t phi)
{
    fPhi = phi / kRadianDeg;
}



void TRadCor::SetMaxMx(Double_t maxMx2)
{
    fMaxMx2 = maxMx2;
}



void TRadCor::Setup(void)
{
    // Calculate the missing mass

    Double_t nu;
    Double_t Sx;

    nu = - fY / ( 2.0 * kMassProton * fX);
    Sx = 2.0 * kMassProton * nu;

    fMx2 = std::pow(kMassProton,2) + Sx * (1.0 - fZ) + fPt;
}



Double_t TRadCor::GetRCFactor(void)
{
    // Get the radiative correction factor. You must set the parameters before
    // using this method.

    return 0.;
}



Double_t TRadCor::GetRCFactor(Double_t Ebeam, Double_t x, Double_t q2, Double_t z,
                              Double_t pt, Double_t phi, Double_t maxMx2)
{
    // Get the radiative correction factor for the given parameters.
    //
    // The Ebeam parameter is the energy of the beam, the x, q2, z, pt and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    SetParameters(Ebeam,x,q2,z,pt,phi,maxMx2);
    return GetRCFactor();
}