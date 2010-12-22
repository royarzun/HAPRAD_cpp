#include "TLorentzInvariants.h"
#include "TRadCor.h"
#include "TGlobalConfig.h"
#include "TKinematicalVariables.h"
#include "THadronKinematics.h"
#include "THapradException.h"

#include "haprad_constants.h"

#include <iostream>


TLorentzInvariants::TLorentzInvariants(const TRadCor* rc)
 : fS(0), fX(0), fSx(0), fSp(0), fQ2(0)
{
    fConfig = rc->GetConfig();
    fKin    = rc->GetKinematicalVariables();
    fHadKin = rc->GetHadronKinematics();
}



TLorentzInvariants::~TLorentzInvariants()
{
    // Do nothing
}



void TLorentzInvariants::Evaluate(void)
{
    using namespace TMath;

    Double_t M = kMassProton;
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

    fS = 2. * M * fKin->E();

    Double_t y;
    if (fKin->Y() >= 0.) {
        y = - fKin->Y();
        fQ2 = fS * fKin->X() * fKin->Y();
    } else {
        fQ2 = - fKin->Y();
        y = fQ2 / (fS * fKin->X());
    }

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "S      " << std::setw(20) << std::setprecision(10) << fS  << std::endl;
    std::cout << "y      " << std::setw(20) << std::setprecision(10) << y   << std::endl;
    std::cout << "Q^2    " << std::setw(20) << std::setprecision(10) << fQ2 << std::endl;
#endif

    Double_t y_max = 1. / (1. + M * M * fKin->X() / fS);
    Double_t y_min = (kMassC2 - M * M) / (fS * (1. - fKin->X()));

    if (y > y_max || y < y_min || fKin->X() > 1. || fKin->X() < 0.)
        throw TKinematicException();

    fX  = fS * (1. - y);
    fSx = fS - fX;
    fSp = fS + fX;
    fW2 = fS - fX - fQ2 + M * M;

    fLambdaS = fS * fS - 4. * m * m * M * M;
    fLambdaX = fX * fX - 4. * m * m * M * M;
    fLambdaM = fQ2 * fQ2 + 4. * m * m * fQ2;
    fLambdaQ = Power(fSx, 2) + 4. * M * M * fQ2;

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "X      " << std::setw(20) << std::setprecision(10) << fX       << std::endl;
    std::cout << "S_x    " << std::setw(20) << std::setprecision(10) << fSx      << std::endl;
    std::cout << "S_p    " << std::setw(20) << std::setprecision(10) << fSp      << std::endl;
    std::cout << "W2     " << std::setw(20) << std::setprecision(10) << fW2      << std::endl;
    std::cout << "l_s    " << std::setw(20) << std::setprecision(10) << fLambdaS << std::endl;
    std::cout << "l_x    " << std::setw(20) << std::setprecision(10) << fLambdaX << std::endl;
    std::cout << "l_m    " << std::setw(20) << std::setprecision(10) << fLambdaM << std::endl;
    std::cout << "l_q    " << std::setw(20) << std::setprecision(10) << fLambdaQ << std::endl;
#endif

    if (fLambdaS < 0.) std::cout << " Conkin: lambda_s < 0 " << std::endl;
    if (fLambdaX < 0.) std::cout << " Conkin: lambda_x < 0 " << std::endl;
    if (fLambdaQ < 0.) std::cout << " Conkin: lambda_q < 0 " << std::endl;
    if (fLambdaM < 0.) std::cout << " Conkin: lambda_m < 0 " << std::endl;

    fSqrtLs = Sqrt(Max(0.,fLambdaS));
    fSqrtLx = Sqrt(Max(0.,fLambdaX));
    fSqrtLq = Sqrt(Max(0.,fLambdaQ));
}



void TLorentzInvariants::EvaluateV12(void)
{
    using namespace TMath;

    Double_t M = kMassProton;
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

    Double_t costs, costx, sints, sintx;
    Double_t lambda;

    costs = (fS * (fS - fX) + 2. * M * M * fQ2) / fSqrtLs / fSqrtLq;
    costx = (fX * (fS - fX) - 2. * M * M * fQ2) / fSqrtLx / fSqrtLq;

    lambda = fS * fX * fQ2 - M * M * fQ2 * fQ2 - m * m * fLambdaQ;

    if (lambda > 0) {
        sints = 2. * M * Sqrt(lambda) / fSqrtLs / fSqrtLq;
        sintx = 2. * M * Sqrt(lambda) / fSqrtLx / fSqrtLq;
    } else {
        std::cout << "sphi: sints = NaN " << lambda << std::endl;
        std::cout << "sphi: sintx = NaN " << lambda << std::endl;
        sints = 0.;
        sintx = 0.;
    }

    Double_t v1, v2;
    v1 = costs * fHadKin->Pl() + sints * fHadKin->Pt() * Cos(fKin->PhiH());
    v2 = costx * fHadKin->Pl() + sintx * fHadKin->Pt() * Cos(fKin->PhiH());

    fV1 = (fS * fHadKin->Eh() - fSqrtLs * v1) / M;
    fV2 = (fX * fHadKin->Eh() - fSqrtLx * v2) / M;

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "V1     " << std::setw(20) << std::setprecision(10) << V_1  << std::endl;
    std::cout << "V2     " << std::setw(20) << std::setprecision(10) << V_2  << std::endl;
#endif
}
