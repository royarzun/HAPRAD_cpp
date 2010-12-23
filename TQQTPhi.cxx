#include "TQQTPhi.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "TRV2LN.h"
#include "haprad_constants.h"
#include "Math/GaussIntegrator.h"


TQQTPhi::TQQTPhi(const TRadCor* rc)
{
    fRC  = rc;
    fInv = rc->GetLorentzInvariants();

    fTauMax = (fInv->Sx() + fInv->SqrtLq()) / (2. * kMassProton * kMassProton);
    fTauMin = - fInv->Q2() / (kMassProton * kMassProton) / fTauMax;

    double tau_1 = - fInv->Q2() / fInv->S();
    double tau_2 =   fInv->Q2() / fInv->X();

    fTauArray[0] = fTauMin;
    fTauArray[1] = tau_1 - 0.15 * (tau_1 - fTauMin);
    fTauArray[2] = tau_1 + 0.15 * (tau_2 - tau_1);
    fTauArray[3] = tau_2 - 0.15 * (tau_2 - tau_1);
    fTauArray[4] = tau_2 + 0.15 * (fTauMax - tau_2);
    fTauArray[5] = fTauMax;
}



TQQTPhi::~TQQTPhi()
{

}



ROOT::Math::IBaseFunctionOneDim* TQQTPhi::Clone() const
{
    return 0;
}



double TQQTPhi::DoEval(double phi) const
{
    ROOT::Math::GaussIntegrator ig;

    TRV2LN rv2ln(fRC, phi);
    ig.SetFunction(rv2ln,false);
    ig.SetRelTolerance(0.01);

    double res = 0;
    Double_t ep = TMath::Power(1, -12);

    for (int i = 0; i < 6; i++) {
        double re = ig.Integral(TMath::Log(fKin->X() + fTauArray[i]) + ep,
                                TMath::Log(fKin->X() + fTauArray[i]) + ep);
        res = res + re;
    }

    return res;
}
