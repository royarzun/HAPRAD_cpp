#include "THadronKinematics.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "THapradException.h"
#include "TLorentzInvariants.h"
#include "haprad_constants.h"

#include <iostream>
#ifdef DEBUG
#include <iomanip>
#endif


THadronKinematics::THadronKinematics(const TRadCor* rc)
 : fEh(0), fPl(0), fPt(0), fNu(0), fPx2(0), fPh(0)
{
    fConfig = rc->GetConfig();
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
}



THadronKinematics::~THadronKinematics()
{
    // Do nothing
}



void THadronKinematics::Evaluate(void)
{
    using namespace TMath;

    Double_t M   = kMassProton;
    Double_t m_h = kMassDetectedHadron;

    fNu = fInv->Sx() / (2. * M);
    fEh = fNu * fKin->Z();
#ifdef DEBUG
    std::cout << "nu     " << std::setw(20) << std::setprecision(10) << fNu << std::endl;
    std::cout << "Eh     " << std::setw(20) << std::setprecision(10) << fEh << std::endl;
#endif

    if (fEh < m_h) throw TKinematicException();

    fPh = Sqrt(fEh * fEh - m_h * m_h);
#ifdef DEBUG
    std::cout << "p_h    " << std::setw(20) << std::setprecision(10) << fPh  << std::endl;
#endif

    fSqNuQ = Sqrt(fNu * fNu + fInv->Q2());

    if (fKin->T() >= 0.) {
        fPt = fKin->T();

        if (fPh < fPt) throw TKinematicException();

        fPl = Sqrt(fPh * fPh - fPt * fPt);
    } else {
        fPl = (fKin->T() + fInv->Q2() - m_h * m_h + 2. * fNu * fEh) / 2. / fSqNuQ;

        if (fPh < Abs(fPl)) {
            Double_t eps1, eps2, eps3, eps4, eps5, sum;

            eps1 = fKin->T() * kEpsMachine / fSqNuQ;
            eps2 = 2. * m_h * m_h  * kEpsMachine / fSqNuQ;
            eps3 = 2. * fNu * fEh * kEpsMachine / fSqNuQ;
            eps4 = fKin->T() + fInv->Q2() - m_h * m_h + 2. * fNu * fEh;
            eps5 = eps4 / fSqNuQ * kEpsMachine;

            sum = eps1 * eps1 + eps2 * eps2 + 2. * eps3 * eps3 + eps5 * eps5;

            Double_t epspl  = Sqrt(sum) / 2.;
            Double_t calEps = fPh - Abs(fPl);
            if (Abs(calEps) > epspl) {
                throw TKinematicException();
            } else {
               std::cout << "Zero p_t! " << fPl
                         << "\t"         << calEps
                         << "\t"         << epspl << std::endl;
               fPl = Sign(1., fPl) * fPh;
            }
        }

        fPt = Sqrt(fPh * fPh - fPl * fPl);
    }
}



void THadronKinematics::EvaluatePx2(void)
{
    using namespace TMath;

    Double_t M   = kMassProton;
    Double_t m_h = kMassDetectedHadron;

    Double_t px2_max = fInv->W2() - m_h * m_h;
    fPx2 = M * M + fInv->Sx() * (1. - fKin->Z()) + fKin->T();

    if (fPx2 < kMassC2 || fPx2 > px2_max) throw TKinematicException();
}
