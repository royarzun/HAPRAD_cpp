#include "TStructFunctionArray.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "THapradUtils.h"
#include "haprad_constants.h"

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
}



TStructFunctionArray::TStructFunctionArray(const TRadCor* rc, Int_t n)
  : TArrayD(n)
{
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();
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

    Double_t tld_Q2  = fInv->Q2() + R * tau;
    Double_t tld_t   = fKin->T() - R * (tau - mu);
    Double_t tld_lq  = (fInv->Sx() - R) * (fInv->Sx() - R) + 4. * M * M * tld_Q2;
    Double_t tld_px2 = fHadKin->Px2() - R * (1. + tau - mu);
    Double_t phq     = (m_h * m_h - tld_Q2 - tld_t) / 2.;

    Double_t tld_nu  = (fInv->Sx() - R) / (2 * M);
    Double_t aks     = tld_Q2 / (2 * M) / tld_nu;
    Double_t zh      = fHadKin->Eh() / tld_nu;
    Double_t tld_sq  = TMath::Sqrt(tld_Q2 + TMath::Power(tld_nu, 2));
    Double_t tld_p_l = (fHadKin->Eh() * tld_nu - phq) / tld_sq;
    Double_t tld_pt2 = TMath::Power(fHadKin->Ph(), 2) - TMath::Power(tld_p_l, 2);


    const Double_t& e = kEpsMachine;

    Double_t eps_nu = Sqrt(
                Power(fInv->Sx() * e / (2 * M), 2) +
                Power(R * e / (2 * M), 2) +
                Power((fInv->Sx() - R) / (2 * M) * e, 2));

    Double_t epst = Sqrt(
                Power((fKin->T() * e), 2) +
                Power(((tau - mu) * R * e), 2) +
                Power((R * tau * e), 2) +
                Power((R * mu * e), 2));

    Double_t epsphq = Sqrt(
                Power((2 * m_h * m_h * e), 2) +
                Power((tld_Q2 * e), 2) +
                Power(epst, 2)
                ) / 2;

    Double_t epsq = 1. / 2. / Sqrt(tld_Q2 + tld_nu * tld_nu) *
                Sqrt(Power((tld_Q2 * e), 2) +
                        Power(2 * tld_nu * tld_nu * eps_nu, 2));

    Double_t epspl = Sqrt(
                Power(fHadKin->Eh() * e * tld_nu / tld_sq, 2) +
                Power(fHadKin->Eh() * eps_nu / tld_sq, 2) +
                Power(epsphq / tld_sq, 2) +
                    Power((fHadKin->Eh() * tld_nu - phq) /
                        Power(tld_sq, 2) * epsq, 2));

    Double_t epspt2 = 2. * Sqrt(
                Power(fHadKin->Ph(), 4) * Power(e, 2) +
                Power(tld_p_l, 4) * Power(epspl, 2));

    if (tld_pt2 >= 0 || (tld_pt2 * tld_pt2 - epspt2 * epspt2) <= 0) {
        Double_t a = fInv->S() / (2 * M) * (fInv->S() / (2 * M) - fHadKin->Nu()) * tld_Q2 / fInv->Q2();
        Double_t tlde = (tld_nu + TMath::Sqrt(TMath::Power(tld_nu, 2) + 4 * a)) / 2;
        Double_t tldy = tld_nu / tlde;

        Double_t H1z, H2z, H3z, H4z;

        HapradUtils::SemiInclusiveModel(tld_Q2, aks, tldy, zh, tld_pt2, tld_px2,
                                        tld_p_l, H1z, H2z, H3z, H4z);

        // Including kinematic shift due to photon emission
        Double_t tld_sqrt_lq = TMath::Sqrt(tld_lq);

        Double_t aa = (fInv->Sx() - R) * (zh - 2 * M * tld_p_l / tld_sqrt_lq) / 2. / (M * M);

        Double_t h1 = zh * (2. * tld_Q2 * aks * H1z - tld_pt2 * H4z) /
                    M / tld_Q2 / aks;

        Double_t h2 = 2 * zh * (2 * aks * (aks * H2z + aa * H3z) * tld_lq +
                H4z * (aa * aa * tld_lq - 2 * tld_pt2 * tld_Q2)) /
                        M / aks / tld_Q2 / tld_lq;

        Double_t h3 = 2 * H4z * zh / M / tld_Q2 / aks;

        Double_t h4 = -2 * zh * (aks * H3z + aa * H4z) / M / tld_Q2 / aks;

        fArray[0] = h1;
        fArray[1] = h2;
        fArray[2] = h3;
        fArray[3] = h4;
    }
}
