#include "TQQTPhi.h"
#include "TRadCor.h"
#include "TGlobalConfig.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "TRV2LN.h"
#include "haprad_constants.h"
#include "Math/GaussLegendreIntegrator.h"
#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif


TQQTPhi::TQQTPhi(const TRadCor* rc)
{
    fRC  = rc;
    fConfig = rc->GetConfig();
    fInv = rc->GetLorentzInvariants();
    fKin = rc->GetKinematicalVariables();

    fTauMax = (fInv->Sx() + fInv->SqrtLq()) / (2. * kMassProton * kMassProton);
    fTauMin = - fInv->Q2() / (kMassProton * kMassProton) / fTauMax;
#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "  tau_max  " << std::setw(20) << std::setprecision(10)
              << fTauMax << std::endl;
    std::cout << "  tau_min  " << std::setw(20) << std::setprecision(10)
              << fTauMin << std::endl;
#endif

    double tau_1 = - fInv->Q2() / fInv->S();
    double tau_2 =   fInv->Q2() / fInv->X();
#ifdef DEBUG
    std::cout << "  tau_1  " << std::setw(20) << std::setprecision(10)
              << tau_1 << std::endl;
    std::cout << "  tau_2  " << std::setw(20) << std::setprecision(10)
              << tau_2 << std::endl;
#endif

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
#ifdef DEBUG
    std::cout << "  QQTPhi(" << phi << ")" << std::endl << std::endl;
#endif
    ROOT::Math::GaussLegendreIntegrator ig;

    TRV2LN rv2ln(fRC, phi);
    ig.SetFunction(rv2ln,false);
    ig.SetNumberPoints(100);
    ig.SetRelTolerance(fConfig->EpsTau());

    double res = 0;
    Double_t ep = TMath::Power(10, -12);

    for (int i = 0; i < 5; i++) {
#ifdef DEBUG
        std::cout.setf(std::ios::fixed);
        std::cout << "  i: " << i << std::endl;
        std::cout << "   lim1 "
                  << std::setw(20) << std::setprecision(10)
                  << TMath::Log(fKin->X() + fTauArray[i]) + ep << std::endl;
        std::cout << "   lim2 "
                  << std::setw(20) << std::setprecision(10)
                  << TMath::Log(fKin->X() + fTauArray[i+1]) + ep << std::endl;
#endif
        double re = ig.Integral(TMath::Log(fKin->X() + fTauArray[i]) + ep,
                                TMath::Log(fKin->X() + fTauArray[i+1]) + ep);
#ifdef DEBUG
        std::cout << "   re   " << std::setw(20) << std::setprecision(10)
                  << re << std::endl;
#endif
        res = res + re;
    }

    return res;
}
