#include "TDelta.h"
#include "TGlobalConfig.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "THapradUtils.h"
#include "haprad_constants.h"


TDelta::TDelta()
  : fVR(0), fInf(0), fVac(0)
{
    // Do nothing
}



TDelta::~TDelta()
{
    // Do nothing
}



void TDelta::Evaluate(const TLorentzInvariants& inv,
                      const THadronKinematics& hkin)
{
    using namespace TMath;

    Double_t m;

    switch(fConfig->PolarizationType()) {
        case 1:
            m = kMassElectron;
            break;
        case 2:
            m = kMassMuon;
            break;
        default:
            m = kMassElectron;
    }

    Double_t S_ = inv.S() - inv.Q2() - inv.V1();
    Double_t X_ = inv.X() + inv.Q2() - inv.V2();

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "S'     " << std::setw(20) << std::setprecision(10) << S_  << std::endl;
    std::cout << "X'     " << std::setw(20) << std::setprecision(10) << X_  << std::endl;
#endif

    Double_t l_m  = Log(inv.Q2() / (m * m));
    Double_t Li_2 = HapradUtils::fspen(1. - hkin.Px2() * inv.Q2() / (S_ * X_));

    fVR  = 1.5 * l_m - 2. - 0.5 * Power(Log(X_ / S_),2) + Li_2 - kPi * kPi / 6.;
    fInf = (l_m - 1.) * Log(Power(hkin.Px2() - kMassC2, 2) / S_ / X_);
    fVac = VacPol(inv.Q2());

    fVR  *= kAlpha / kPi;
    fInf *= kAlpha / kPi;
    fVac *= kAlpha / kPi;
}



Double_t TDelta::VacPol(const Double_t Q2)
{
    Double_t leptonMass[3] = { 0.26110  * TMath::Power(10,-6),
                               0.111637 * TMath::Power(10,-1),
                               3.18301 };

    Double_t suml = 0;
    for (Int_t i = 0; i < 3; ++i) {
        Double_t a2    = 2. * leptonMass[i];
        Double_t sqlmi = TMath::Sqrt(Q2 * Q2 + 2. * a2 * Q2);
        Double_t allmi = TMath::Log((sqlmi + Q2) / (sqlmi - Q2)) / sqlmi;

        suml = suml + 2. * (Q2 + a2) * allmi / 3. - 10. / 9. +
                    4. * a2 * (1. - a2 * allmi) / 3. / Q2;
    }

    Double_t a, b, c;

    if (Q2 < 1) {
        a = -1.345 * TMath::Power(10,-9);
        b = -2.302 * TMath::Power(10,-3);
        c = 4.091;
    } else if (Q2 < 64) {
        a = -1.512 * TMath::Power(10,-3);
        b = -2.822 * TMath::Power(10,-3);
        c =  1.218;
    } else {
        a = -1.1344 * TMath::Power(10,-3);
        b = -3.0680 * TMath::Power(10,-3);
        c =  9.9992 * TMath::Power(10,-1);
    }

    Double_t sumh;
    sumh = - (a + b * TMath::Log(1. + c * Q2)) * 2. * kPi / kAlpha;

#ifdef DEBUG
    std::cout << std::endl;
    std::cout.setf(std::ios::fixed);
    std::cout << "suml   " << std::setw(20) << std::setprecision(10) << suml  << std::endl;
    std::cout << "sumh   " << std::setw(20) << std::setprecision(10) << sumh  << std::endl;
    std::cout << std::endl;
#endif

    return suml + sumh;
}
