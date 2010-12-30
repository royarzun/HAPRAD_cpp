#include "TSffun.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "TMath.h"
#include "TExclusiveModel.h"
#include "haprad_constants.h"
#include <iostream>


TSffun::TSffun(const TRadCor* rc)
    : TArrayD(4)
{
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();

}

TSffun::~TSffun()
{
    // Do Nothing, default destructor
}

void TSffun::Evaluate(Double_t Q2, Double_t w2, Double_t t)
{
    using namespace TMath;

    Double_t st, sl, stt, slt, sltp, sfm10;
    Double_t M_p = kMassProton;
    Double_t m = kMassDetectedHadron;
    Double_t M_u = kMassUndetectedHadron;
    Double_t SqrtW2 = Sqrt(w2);
    Double_t Sx_t = fInv->Sx() + t + M_p * M_p - M_u*M_u;
    Double_t tq = t + Q2 - m * m;
    Double_t sffun_cmp = Power((w2 - M_u*M_u - m * m), 2) - 4 * M_u*M_u * m * m;

    if (sffun_cmp < 0)
        std::cout << "sffun: SqrtLw=NaN " << sffun_cmp << std::endl;
    

    Double_t SqrtLw = Sqrt(Max(0., sffun_cmp));
    Double_t ssffun_cmp = Q2 * Power(Sx_t, 2) - Sx_t * fInv->Sx() * tq - M_p * M_p * Power(tq, 2) - m * m * fInv->LambdaQ();

    if (ssffun_cmp < 0)
        std::cout << "ssffun: qll=NaN " << ssffun_cmp << std::endl;
    

    Double_t SqrtLl = Sqrt(Max(0., ssffun_cmp));
    Double_t cspion = (2 * tq * w2 + (fInv->Sx() - 2 * Q2) * (w2 + m * m - M_u*M_u)) / SqrtLw / Sqrt(fInv->LambdaQ());

//  Exclusive peak model (cross sections sigma_L,T,LT... from MAID2003)

    ExclusiveModel(Q2, SqrtW2, cspion, st, sl, stt, slt, sltp);

    Double_t sfm20;
    Double_t sfm2tl;
    Double_t sfm4tl;
    Double_t sfm5tl;
    Double_t sfm4tt;
    Double_t sfm3tt;
    Double_t sfm2tt;
    Double_t Coetr;
// Structure functions

    if (fInv->LambdaQ() > 0 && SqrtLl > 0 && SqrtLw > 0) {
        sfm10 = st - stt;
        sfm20 = 4. * (st + sl) * Q2 / fInv->LambdaQ();
        sfm2tl = 2. * slt * Sqrt(Q2) * (-fInv->Sx() * tq + 2. * Q2 * Sx_t) / (fInv->LambdaQ() * SqrtLl);
        sfm4tl = -slt * Sqrt(Q2) / SqrtLl;
        sfm4tt = -2. * stt * (-fInv->Sx() * tq + 2. * Q2 * Sx_t) / Power(SqrtLl, 2.);
        sfm3tt = 2. * stt * fInv->LambdaQ() / Power(SqrtLl, 2.);
        sfm2tt = 2. * stt * (Power((-fInv->Sx() * tq + 2. * Q2 * Sx_t), 2.) - 2. * Q2 * Power(SqrtLl, 2.)) / (fInv->LambdaQ() * Power(SqrtLl, 2.));
        sfm5tl = -sltp * Sqrt(Q2) / SqrtLl;

        Coetr = 16. * kPi * (w2 - M_p * M_p) * w2 / (kAlpha * SqrtLw) / kBarn * 1000.;
        fArray[0] = Coetr * sfm10;
        fArray[1] = Coetr * (sfm20 + sfm2tl + sfm2tt);
        fArray[2] = Coetr * sfm3tt;
        fArray[3] = Coetr * (sfm4tl + sfm4tt);

    } else {
        fArray[0] = 0;
        fArray[1] = 0;
        fArray[2] = 0;
        fArray[3] = 0;
    }

    if (fArray[2] != fArray[2])
        std::cout << "sffun: "      << Coetr         << st
                  << sl             << stt           << slt
                  << sltp           << Q2            << fInv->LambdaQ()
                  << fInv->Sx()     << SqrtW2        << cspion
                  << Sx_t           << SqrtLl          << std::endl;

}


