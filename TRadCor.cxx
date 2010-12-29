#include "TRadCor.h"
#include "THapradException.h"
#include "TDelta.h"
#include "TQQTPhi.h"
#include "TBorn.h"
#include "haprad_constants.h"
#include "Math/GaussLegendreIntegrator.h"
#include <iostream>
#ifdef DEBUG
#include <iomanip>
#endif


TRadCor::TRadCor()
: fInv(this), fHadKin(this),
  maxMx2(0.), Mx2(0.),
  sigma_born(0.), sig_obs(0.), tail(0.),
  M(kMassProton), m(kMassElectron), m_h(kMassDetectedHadron)
{
    // Default constructor
}



TRadCor::TRadCor(Double_t E, Double_t x, Double_t Q2, Double_t z,
                 Double_t p_t, Double_t phi, Double_t maxMx2)
: fInv(this), fHadKin(this),
  maxMx2(0.), Mx2(0.),
  sigma_born(0.), sig_obs(0.), tail(0.),
  M(kMassProton), m(kMassElectron), m_h(kMassDetectedHadron)
{
    // Normal constructor for a radiative correction object
    //
    // The E parameter is the energy of the beam, the x, Q2, z, p_t and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount
    // of missing mass.

    fKin.SetAll(x,-Q2,z,p_t,phi/kRadianDeg,E);
    Setup();
}



TRadCor::~TRadCor()
{
    // Default destructor
}



void TRadCor::SetParameters(Double_t E, Double_t x, Double_t Q2, Double_t z,
                            Double_t p_t, Double_t phi, Double_t maxMx2)
{
    // Set the values for the variables used to calculate the radiative
    // correction.
    //
    // The E parameter is the energy of the beam, the x, Q2, z, p_t and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount
    // of missing mass.

    this->maxMx2 = maxMx2;

    fKin.SetAll(x,-Q2,z,p_t,phi/kRadianDeg,E);
    Setup();
}



void TRadCor::SetEbeam(Double_t E)
{
    fKin.SetE(E);
}



void TRadCor::SetX(Double_t x)
{
    fKin.SetX(x);
    Setup();
}



void TRadCor::SetQ2(Double_t Q2)
{
    fKin.SetY(-Q2);
    Setup();
}



void TRadCor::SetZ(Double_t z)
{
    fKin.SetZ(z);
    Setup();
}



void TRadCor::SetPt(Double_t p_t)
{
    fKin.SetT(p_t);
    Setup();
}



void TRadCor::SetPhi(Double_t phi)
{
    fKin.SetPhiH(phi/kRadianDeg);
}



void TRadCor::SetMaxMx(Double_t maxMx2)
{
    this->maxMx2 = maxMx2;
}



void TRadCor::RegisteredLepton(Int_t type)
{
    // Set the registere lepton: 1 -- electron; 2 -- muon.
    //
    // Default is 1.

    fConfig.SetLepton(type);
}



void TRadCor::IntegratePhiRad(Int_t type)
{
    // Set whether to integrate over phi_{rad} (1) or to approximate (0).
    //
    // Default is 0;

    fConfig.SetIntegrationPhiRad(type);
}



void TRadCor::IntegratePhiHad(Int_t type)
{
    // Set whether to integrate over phi_{had} (1) or not (0).
    //
    // Default is 0;

    fConfig.SetIntegrationPhiHad(type);
}



void TRadCor::SetPolarization(Int_t type)
{
    // Set the type of the polarization:
    //
    //     1 -- long; 2 -- tran; 0 -- unpol.
    //
    // Default is 0;

    fConfig.SetPolarization(type);
}



void TRadCor::Setup(void)
{
    // Calculate the missing mass squared

    Double_t S_x = - fKin.Y() / fKin.X();
    Mx2 = M * M + S_x * (1.0 - fKin.Z()) + fKin.T();
}



Double_t TRadCor::GetRCFactor(void)
{
    // Get the radiative correction factor. You must set the parameters before
    // using this method.

    if (Mx2 > maxMx2) {
        Haprad();
        rc = sig_obs / sigma_born;
    } else {
        rc = 0;
    }
    return rc;
}



Double_t TRadCor::GetRCFactor(Double_t E, Double_t x, Double_t Q2, Double_t z,
                              Double_t p_t, Double_t phi, Double_t maxMx2)
{
    // Get the radiative correction factor for the given parameters.
    //
    // The E parameter is the energy of the beam, the x, Q2, z, p_t and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    SetParameters(E,x,Q2,z,p_t,phi,maxMx2);
    return GetRCFactor();
}



void TRadCor::Haprad(void)
{
    pl = 0.;

    try {
        fInv.Evaluate();
        if (fKin.Y() < 0.) {
            Double_t y = fInv.Q2() / (fInv.S() * fKin.X());
            fKin.SetY(y);
        }
    } catch (TKinematicException& wrongKin) {
        std::cout << wrongKin.what() << std::endl;
        return;
    }

    try {
        fHadKin.Evaluate();
    } catch (TKinematicException& wrongKin) {
        std::cout << wrongKin.what() << std::endl;
        return;
    }


    N = kPi * TMath::Power(kAlpha,2) * fKin.Y() * fInv.Sx() * M / 2. / fInv.SqrtLq() * kBarn;
#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "N      " << std::setw(20)
                           << std::setprecision(10) << N << std::endl;
#endif

    if (fKin.T() >= 0.) {
        if (fHadKin.Ph() > fHadKin.Pt()) {
            N = N * fInv.SqrtLq() / 2. / M / fHadKin.Pl();
        } else {
            N = 0.;
        }

        Double_t t = m_h * m_h - fInv.Q2() + 2. * (fHadKin.SqNuQ() * fHadKin.Pl() - fHadKin.Nu() * fHadKin.Eh());
        fKin.SetT(t);

        std::cout << "p_l: " << fHadKin.Pl() << "\t" << t << "\t"
                  << fHadKin.Pl() - t + fInv.Q2() - m_h * m_h + 2. * fHadKin.Nu() * fHadKin.Eh() / 2. / fHadKin.SqNuQ()
                  << std::endl;
    }
#ifdef DEBUG
    std::cout << "Pt     " << std::setw(20) << std::setprecision(10)
                           << fHadKin.Pt() << std::endl;
    std::cout << "Pl     " << std::setw(20) << std::setprecision(10)
                           << fHadKin.Pl() << std::endl;
    std::cout << "tdif   " << std::setw(20) << std::setprecision(10)
                           << fKin.T() << std::endl;
#endif

    fInv.EvaluateV12();
    fHadKin.EvaluatePx2();

    t_min = m_h * m_h - fInv.Q2() + 2. * (fHadKin.SqNuQ() * fHadKin.Ph() - fHadKin.Nu() * fHadKin.Eh());

    Double_t tdmax;
    tdmax = m_h * m_h - fInv.Q2() + 2. * (- fHadKin.SqNuQ() * fHadKin.Ph() - fHadKin.Nu() * fHadKin.Eh());

    try {
        if ((fKin.T() - t_min) > kEpsMachine || fKin.T() < tdmax) {
            throw TKinematicException();
        }
    } catch (TKinematicException& wrongKin) {
        std::cout << wrongKin.what()
                << " t     = " << fKin.T() << std::endl
                << " t_max = " << tdmax << std::endl
                << " t_min = " << t_min << " - " << fKin.T() - t_min << std::endl;
        return;
    }

    SPhiH();
}



void TRadCor::SPhiH(void)
{
    TDelta fDeltas(this);
    fDeltas.Evaluate();

    Double_t sibt;
    Int_t it_end = 3;

    if (fConfig.PolarizationType() == 0) {
        it_end = 1;
    }

    std::cout << "********** ita: " << ita
            << " *********" << std::endl;
    TBorn fBornin(this);
    sigma_born = fBornin.GetValue(N);
    BorninTest(sibt);
    std::cout << "sib1" << sigma_born << std::endl;
    std::cout << "sibt" << sibt << std::endl;
    if (sigma_born == 0.0) {
        tai[1] = 0.;
    }
    std::cout << "tai[" << 1
            << "]\t"  << tai[1] << std::endl;

    qqt(tai[1]);

    Double_t extai1 = TMath::Exp(fDeltas.Inf());
    sig_obs = sigma_born * extai1 * (1. + fDeltas.VR() + fDeltas.Vac()) +
                                                            tai[1] + tai[2];
}



void TRadCor::BorninTest(Double_t& sigma_born)
{

}



void TRadCor::qqt(Double_t& tai)
{
    /*
    if (ita == 1) {
    */

    TQQTPhi qphi(this);
    if (fConfig.IntegratePhiRad() == 1) {
        ROOT::Math::GaussLegendreIntegrator ig;
        ig.SetFunction(qphi);
        ig.SetNumberPoints(150);
        ig.SetRelTolerance(fConfig.EpsPhiR());
        tai = ig.Integral(0, TMath::TwoPi());
        tai = N * kAlpha * tai / (kPi * kPi) / 4. / fInv.SqrtLq();
    } else if (fConfig.IntegratePhiRad() == 0) {
        tai = N * kAlpha / kPi * qphi(0.) / 2 / fInv.SqrtLq();
    }
    /*
    } else {
        Double_t tau_1, tau_2;
        Double_t tau[6];
        Double_t phi[4];

        tau_1 = - Q2 / S;
        tau_2 =   Q2 / X;

        phi[0] = 0.;
        phi[1] = 0.01 * kPi;
        phi[2] = 2. * kPi - 0.01 * kPi;
        phi[3] = 2. * kPi;

        tau[0] = tau_min;
        tau[1] = tau_1 - 0.15 * (tau_1 - tau_min);
        tau[2] = tau_1 + 0.15 * (tau_2 - tau_1);
        tau[3] = tau_2 - 0.15 * (tau_2 - tau_1);
        tau[4] = tau_2 + 0.15 * (tau_max - tau_2);
        tau[5] = tau_max;

        Int_t id = 1;
        Double_t rere = 0.;

        for (Int_t iph = 0; iph < 3; iph++) {
            for (Int_t ico = 0; ico < 5; ico++) {
                Double_t am[2], bm[2];

                am[0] = tau[ico];
                bm[0] = tau[ico+1];
                am[1] = phi[iph];
                bm[1] = phi[iph+1];

                if (am[1] > bm[1])
                    std::cout << am[1] << " < " << bm[1] << std::endl;
                if (am[0] > bm[0])
                    std::cout << am[0] << " < " << bm[1] << std::endl;

                Double_t otr;
//                Double_t ot = TMath::Power(10, -3);
                Int_t mir = 10000;
//                Int_t ma = 10 * mir;
//                Double_t wrk[500];
                Double_t re;
//                Double_t re2;

                // Integration using d01fce
//              d01fce(2, am, bm, mir, ma, rv2tr, ot, otr, 500, wrk, re, id);
                std::cout.setf(std::ios::fixed);
                std::cout << " tai: "
                          << std::setw(4)  << ico + 1
                          << std::setw(4)  << iph + 1
                          << std::setw(10) << std::setprecision(4) << re
                          << std::setw(10) << std::setprecision(4) << otr
                          << std::setw(10) << std::setprecision(4) << mir
                          << std::setw(10) << std::setprecision(4) << id
                          << std::endl;

                rere = rere + re;
            }
        }

        tai = - kAlpha / (64. * TMath::Power(kPi,5.) * sqrt_lq * M) * N * rere;
*/
}
