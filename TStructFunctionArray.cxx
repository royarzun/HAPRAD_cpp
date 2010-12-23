#include "TStructFunctionArray.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "THapradUtils.h"
#include "haprad_constants.h"

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif

//______________________________________________________________________________
//
// The function calculates deep inelastic (ita = 1), elastic (ita = 2),
// quasielastic (ita = 3) structure functions in kinematical point (tau,R).
//
//     R = S_x - tt
//     tau = (t - Q2) / R
//
// where
//
//     tt = t + mf2 - M^2
//
// mf2 is invarint mass of final hadrons


TStructFunctionArray::TStructFunctionArray(const TRadCor* rc)
  : TArrayD(4)
{
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();
#ifdef DEBUG
    std::cout << "  Structure Function " << std::endl;
#endif
}



TStructFunctionArray::TStructFunctionArray(Int_t n, const TRadCor* rc)
  : TArrayD(n)
{
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();
#ifdef DEBUG
    std::cout << "  Structure Function " << std::endl;
#endif
}



TStructFunctionArray::~TStructFunctionArray()
{
    // Do nothing
}



void TStructFunctionArray::Evaluate(Double_t tau, Double_t mu, Double_t R)
{
    using namespace TMath;

    const Double_t& M = kMassProton;
    const Double_t& m_h = kMassDetectedHadron;

    Double_t tldQ2  = fInv->Q2() + R * tau;
    Double_t tldSx  = fInv->Sx() - R ;

    Double_t tldLq  = tldSx * tldSx + 4. * M * M * tldQ2;
    Double_t tldNu  = tldSx / (2 * M);

    Double_t tldT   = fKin->T() - R * (tau - mu);

    Double_t tldPx2 = fHadKin->Px2() - R * (1. + tau - mu);
    Double_t PhQ    = (m_h * m_h - tldQ2 - tldT) / 2.;

    Double_t tldX   = tldQ2 / (2 * M * tldNu);
    Double_t tldZ   = fHadKin->Eh() / tldNu;

    Double_t tld_sq = Sqrt(tldQ2 + tldNu * tldNu);
    Double_t tldPl  = (fHadKin->Eh() * tldNu - PhQ) / tld_sq;
    Double_t tldPt2 = Power(fHadKin->Ph(), 2) - Power(tldPl, 2);



    Double_t epsNu = kEpsMachine * Sqrt(Power(fHadKin->Nu(), 2) +
                                        Power(fHadKin->Nu() - tldNu, 2) +
                                        Power(tldNu, 2));

    Double_t epsT = kEpsMachine * Sqrt(Power(fKin->T(), 2) +
                                       Power((tau - mu) * R, 2) +
                                       Power(R * tau, 2) +
                                       Power(R * mu, 2));

    Double_t epsPhq = Sqrt(Power(2 * m_h * m_h * kEpsMachine, 2) +
                           Power(tldQ2 * kEpsMachine, 2) +
                           Power(epsT, 2)) / 2;

    Double_t epsQ2 = Sqrt(Power(tldQ2 * kEpsMachine, 2) +
                          Power(2 * tldNu * tldNu * epsNu, 2)) /
                     (2. * Sqrt(tldQ2 + Power(tldNu, 2)));

    Double_t epsPl = Sqrt(Power(fHadKin->Eh() * kEpsMachine * tldNu / tld_sq, 2) +
                          Power(fHadKin->Eh() * epsNu / tld_sq, 2) +
                          Power(epsPhq / tld_sq, 2) +
                          Power((fHadKin->Eh() * tldNu - PhQ) /
                                        Power(tld_sq, 2) * epsQ2, 2));

    Double_t epsPt2 = 2. * Sqrt(Power(fHadKin->Ph(), 4) * Power(kEpsMachine, 2) +
                                Power(tldPl, 4) * Power(epsPl, 2));

    if (tldPt2 >= 0 || (tldPt2 * tldPt2 - epsPt2 * epsPt2) <= 0) {
        Double_t a = fInv->S() / (2 * M) *
                     (fInv->S() / (2 * M) - fHadKin->Nu()) *
                     tldQ2 / fInv->Q2();

        Double_t tldE = 0.5 * (tldNu + Sqrt(Power(tldNu, 2) + 4 * a));
        Double_t tldY = tldNu / tldE;

        Double_t H1z, H2z, H3z, H4z;

        HapradUtils::SemiInclusiveModel(tldQ2, tldX, tldY, tldZ, tldPt2,
                                        tldPx2, tldPl, H1z, H2z, H3z, H4z);

        // Including kinematic shift due to photon emission
        Double_t aa = tldNu * (tldZ - 2 * M * tldPl / Sqrt(tldLq)) / (2. * M * M);

        Double_t h1 = tldZ * (2. * tldQ2 * tldX * H1z - tldPt2 * H4z) /
                    (M * tldQ2 * tldX);

        Double_t h2 = 2 * tldZ * (2 * tldX * (tldX * H2z + aa * H3z) * tldLq +
                H4z * (aa * aa * tldLq - 2 * tldPt2 * tldQ2)) /
                        M / tldX / tldQ2 / tldLq;

        Double_t h3 = 2 * H4z * tldZ / M / tldQ2 / tldX;

        Double_t h4 = -2 * tldZ * (tldX * H3z + aa * H4z) / M / tldQ2 / tldX;

        fArray[0] = h1;
        fArray[1] = h2;
        fArray[2] = h3;
        fArray[3] = h4;
    }
}
