#include "TPODINL.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "TStructFunctionArray.h"
#include "TThetaMatrix.h"
#include "haprad_constants.h"
#include "Math/GaussIntegrator.h"


TPODINL::TPODINL(const TRadCor* rc, double tau, double mu,
                 const TStructFunctionArray& H0, const TThetaMatrix& theta)
  : fTau(tau), fMu(mu), fH0(H0), fTheta(theta)
{
    fRC     = rc;
    fInv    = rc->GetLorentzInvariants();
#ifdef DEBUG
    std::cout << "  PODINL " << std::endl;
#endif
}



TPODINL::~TPODINL()
{

}



ROOT::Math::IBaseFunctionOneDim* TPODINL::Clone() const
{
    return 0;
}



double TPODINL::DoEval(double R) const
{
    double pp, pres;

    TStructFunctionArray H(fRC);
    H.Evaluate(fTau, fMu, R);

    double podinl = 0;
    for (int isf = 0 ; isf < 4; isf++) {
        for (int irr = 0; irr < 2; irr++) {
            pp = H[isf];
            if (irr == 1) pp = pp - fH0[isf] * TMath::Power((1 + R * fTau / fInv->Q2()), 2);

            pres = pp * TMath::Power(R, (irr - 1)) / TMath::Power((fInv->Q2() + R * fTau), 2);
            podinl = podinl - fTheta[isf][irr] * pres;
        }
    }
    return podinl;
}
