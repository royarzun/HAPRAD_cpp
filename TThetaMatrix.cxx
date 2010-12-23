#include "TThetaMatrix.h"
#include "TRadCor.h"
#include "TGlobalConfig.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "THapradUtils.h"
#include "haprad_constants.h"

#include <iostream>

TThetaMatrix::TThetaMatrix(const TRadCor* rc)
  : TMatrixD(4,3)
{
    fConfig = rc->GetConfig();
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();
}



TThetaMatrix::TThetaMatrix(Int_t rows, Int_t cols, const TRadCor* rc)
  : TMatrixD(rows,cols)
{
    fConfig = rc->GetConfig();
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();
}



TThetaMatrix::~TThetaMatrix()
{

}



void TThetaMatrix::Evaluate(Double_t tau, Double_t mu,
                            Int_t ita, Double_t phi_k)
{
    using namespace TMath;

    const Double_t M   = kMassProton;
    const Double_t M2  = kMassProton * kMassProton;
    const Double_t mh2 = kMassDetectedHadron * kMassDetectedHadron;

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
    Double_t m2 = m * m;

    Double_t bb;
    Double_t sqrtmb;
    Double_t z1;
    Double_t z2;
    Double_t bi12;
    Double_t bi1pi2;
    Double_t bis;
    Double_t bir;
    Double_t b1i;
    Double_t b11i;

    // Use simplified names
    const Double_t& S  = fInv->S();
    const Double_t& X  = fInv->X();
    const Double_t& Sx = fInv->Sx();
    const Double_t& Sp = fInv->Sp();
    const Double_t& Q2 = fInv->Q2();
    const Double_t& Lq = fInv->LambdaQ();
    const Double_t& sqrtLq = fInv->SqrtLq();

    Double_t tau2 = tau * tau;

    if (fConfig->IntegratePhiRad() == 0 && ita == 1) {
        Double_t b1;
        Double_t b2;
        Double_t c1;
        Double_t c2;

        b1 = (- Lq * tau - Sp * Sx * tau - 2. * Sp * Q2) / 2.;
        b2 = (- Lq * tau + Sp * Sx * tau + 2. * Sp * Q2) / 2.;

        c1 = Power((S * tau + Q2), 2) - (4. * (M2 * tau2 - Sx * tau - Q2) * m2);
        c2 = Power((S * tau - Q2), 2) - (4. * (M2 * tau2 - Sx * tau - Q2) * m2);

        bb = 1;
        if (c1 < 0) {
            std::cout << "tails: sc1=NaN " << c1 << std::endl;
        }
        Double_t sc1 = TMath::Sqrt(TMath::Max(0., c1));

        if (c2 < 0) {
            std::cout << "tails: sc2=NaN " << c2 << std::endl;
        }
        Double_t sc2 = TMath::Sqrt(TMath::Max(0., c2));

        bi12 = sqrtLq * (Sp * (Sx * tau + 2. * Q2)) / (sc1 * sc2 * (sc1 + sc2));
        bi1pi2 = sqrtLq / sc2 + sqrtLq / sc1;

        bis  =   sqrtLq * (-b1 / sc1 / c1 + b2 / sc2 / c2) * m2;
        bir  =   sqrtLq * (b2 / sc2 / c2 + b1 / sc1 / c1) * m2;
        b1i  = - sqrtLq * b1 / Lq / sqrtLq;
        b11i =   sqrtLq * (3. * b1 * b1 - Lq * c1) / 2. / (Lq * Lq) / sqrtLq;
    } else {
        Double_t tau_max = (Sx + sqrtLq) / (2. * M * M);
        Double_t tau_min = - Q2 / (M * M) / tau_max;
        Double_t Lt = (tau - tau_min) * (tau_max - tau);

        const Double_t Q4 = Q2 * Q2;

        Double_t sqrtmb_comp =  Lt * (S * X * Q2 - Q4 * M2 - m2 * Lq);

        if ( sqrtmb_comp > 0) {
            sqrtmb = TMath::Sqrt(sqrtmb_comp);
        } else {
            std::cout << "tails: sqrtmb=NaN " << sqrtmb_comp << std::endl;
            sqrtmb = 0;
        }

        z1 = (Q2 * Sp + tau * (S * Sx + 2 * M2 * Q2) -
                                2 * M * Cos(phi_k) * sqrtmb) / Lq;
        z2 = (Q2 * Sp + tau * (X * Sx - 2 * M2 * Q2) -
                                2 * M * Cos(phi_k) * sqrtmb) / Lq;
        bb     = 1. / sqrtLq / kPi;
        bi12   = bb / (z1 * z2);
        bi1pi2 = bb / z2 + bb / z1;
        bis    = (bb / Power(z2, 2) + bb / Power(z1, 2)) * m2;
        bir    = (bb / Power(z2, 2) - bb / Power(z1, 2)) * m2;
        b1i    = bb * z1;
        b11i   = bb * Power(z1, 2);
    }

    const Double_t& zh = fKin->T();
    const Double_t& V1 = fInv->V1();
    const Double_t& V2 = fInv->V2();

    Double_t hi2 = bis - (Q2 + 2 * m2) * bi12;
    Double_t vvp = (V1 + V2) / 2.;
    Double_t vvm = (V1 - V2) / 2.;

    (*this)[0][0] = 4. * Q2 * hi2;
    (*this)[0][1] = 4. * tau * hi2;
    (*this)[0][2] = -2. * (bi12 * tau2 + 2. * bb);

    (*this)[1][0] = 2. * hi2 * (S * X - M2 * Q2);
    (*this)[1][1] = 0.5 * (bi1pi2 * Sx * Sp - bi12 * tau * Sp * Sp +
                            2. * bir * Sp + 2. * hi2 * (Sx - 2. * M2 * tau));
    (*this)[1][2] = 0.5 * (bi12 * tau * (2. * M2 * tau - Sx) -
                                            bi1pi2 * Sp + 4. * M2 * bb);

    (*this)[2][0] = 1. * hi2 * (V1 * V2 - mh2 * Q2);
    (*this)[2][1] = - 2. * ((mh2 * tau - mu * vvm) * hi2 - bir * mu * vvp -
                            bi1pi2 * vvm * vvp + bi12 * tau * vvp * vvp);
    (*this)[2][2] = bi12 * tau * (mh2 * tau - mu * vvm) -
                    bi1pi2 * mu * vvp + 2. * mh2 * bb;

    (*this)[3][0] = -1. * (V1 * Sx - 2. * S * vvp + Q2 * Sx * zh) * hi2;
    (*this)[3][1] = - 2. * bi12 * tau * vvp * Sp +
                    bi1pi2 * (Sx * vvm + Sp * vvp) +
                    bir * (mu * Sp + 2. * vvp) +
                    hi2 * (Sx * (mu - 2. * tau * zh) + 2. * vvm);
    (*this)[3][2] = bi12 * tau * (Sx * (tau * zh - mu / 2.) - vvm) -
                    bi1pi2 * (Sp * mu / 2. + vvp) + 2. * Sx * zh * bb;
}
