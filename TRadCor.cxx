#include "TRadCor.h"
#include "THapradException.h"
#include "TDelta.h"
#include "haprad_constants.h"
#include <iostream>
#include <iomanip>


TRadCor::TRadCor()
: fInv(this), fHadKin(this),
  maxMx2(0.), Mx2(0.),
  sigma_born(0.), sig_obs(0.), tail(0.),
  M(kMassProton), m(kMassElectron), m_h(kMassDetectedHadron),
  eps_phir(0.01), eps_tau(0.001), eps_rr(0.001)
{
    // Default constructor
}



TRadCor::TRadCor(Double_t E, Double_t x, Double_t Q2, Double_t z,
                 Double_t p_t, Double_t phi, Double_t maxMx2)
: fInv(this), fHadKin(this),
  maxMx2(0.), Mx2(0.),
  sigma_born(0.), sig_obs(0.), tail(0.),
  M(kMassProton), m(kMassElectron), m_h(kMassDetectedHadron),
  eps_phir(0.01), eps_tau(0.001), eps_rr(0.001)
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
    std::cout << "N      " << std::setw(20) << std::setprecision(10) << N << std::endl;
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
    std::cout << "Pt     " << std::setw(20) << std::setprecision(10) << fHadKin.Pt() << std::endl;
    std::cout << "Pl     " << std::setw(20) << std::setprecision(10) << fHadKin.Pl() << std::endl;
    std::cout << "tdif   " << std::setw(20) << std::setprecision(10) << fKin.T() << std::endl;
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

    for (Int_t i = 1; i <= it_end; ++i) {
        for (ita = 1; ita <= 2; ++ita) {
            std::cout << "********** ita: " << ita
                      << " *********" << std::endl;
            if (ita == 1) {
                sigma_born = Bornin();
                BorninTest(sibt);
                std::cout << "sib1" << sigma_born << std::endl;
                std::cout << "sibt" << sibt << std::endl;
                if (sigma_born == 0.0) {
                    tai[1] = 0.;
                    continue;
                }
            }
            qqt(tai[ita]);
            std::cout << "tai[" << ita
                      << "]\t"  << tai[ita] << std::endl;
        }

        Double_t extai1 = TMath::Exp(fDeltas.Inf());
        sig_obs = sigma_born * extai1 * (1. + fDeltas.VR() + fDeltas.Vac()) +
                                                            tai[1] + tai[2];
    }
}






Double_t TRadCor::Bornin(void)
{
    /*
    Double_t H[4];
    Double_t thetaB[4];

    strf(0.,0.,0,H);

    thetaB[0] = Q2;
    thetaB[1] = (S * X - M * M * Q2) / 2.;
    thetaB[2] = (V_1 * V_2 - m_h * m_h * Q2) / 2.;
    thetaB[3] = (V_2 * S + V_1 * X - z * Q2 * S_x) / 2.;

#ifdef DEBUG
        std::cout << "    BORNIN      " << std::endl;
#endif

    Double_t sum = 0.;
    for (Int_t i = 0; i < 4; ++i) {
        sum = sum + thetaB[i] * H[i];
#ifdef DEBUG
        std::cout.setf(std::ios::fixed);
        std::cout << "    i         " << std::setw(20) << std::setprecision(10) << i         << std::endl;
        std::cout << "    theta^B_i " << std::setw(20) << std::setprecision(10) << thetaB[i] << std::endl;
        std::cout << "    H[i]      " << std::setw(20) << std::setprecision(10) << H[i]      << std::endl;
        std::cout << "    sum       " << std::setw(20) << std::setprecision(10) << sum       << std::endl;
#endif
    }

    return sum * N / Q2 / Q2 * 2.;
*/
    return 0;
}



void TRadCor::BorninTest(Double_t& sigma_born)
{

}



void TRadCor::qqt(Double_t& tai)
{
    /*
    Double_t sqrt_lq = TMath::Sqrt(TMath::Max(0., lambda_q));

    if (ita == 1) {
        if (int_phi_rad == 1) {
            // Integration using Simpson's rule
//            simpsx(0, (2. * pi), 150, epsphir, qqtphi, tai);
            tai = N * kAlpha * tai / (kPi * kPi) / 4. / sqrt_lq;
        } else if (int_phi_rad == 0) {
            tai = N * kAlpha / kPi * qqtphi(0) / 2 / sqrt_lq;
        }
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
    }
*/
}



Double_t TRadCor::qqtphi(Double_t phi)
{
/*
    phi_rad = phi;

    Double_t tau_1, tau_2;
    Double_t tau[6];

    tau_1 = - Q2 / S;
    tau_2 =   Q2 / X;

    tau[0] = tau_min;
    tau[1] = tau_1 - 0.15 * (tau_1 - tau_min);
    tau[2] = tau_1 + 0.15 * (tau_2 - tau_1);
    tau[3] = tau_2 - 0.15 * (tau_2 - tau_1);
    tau[4] = tau_2 + 0.15 * (tau_max - tau_2);
    tau[5] = tau_max;

//    Double_t ep = TMath::Power(1, -12);
    Double_t res = 0;

    for (Int_t i = 0; i < 6; i++) {
        Double_t re;
//        simptx(TMath::Log(xs + tar[i]) + ep, TMath::Log(xs + tar[i+1]) - ep, 100, epstau, rv2ln, re));
        res = res + re;
    }
    return res;
*/
    return 0;
}



Double_t TRadCor::rv2ln(Double_t tauln)
{
    /*
    Double_t tau, mu, factor;
    Double_t tm[4][3];
    tau = TMath::Exp(tauln) - Q2 / S_x;

    Double_t sqrt_lq = TMath::Sqrt(TMath::Max(0., lambda_q));
    Double_t costk, sintk;
    costk = (S_x - M * M * tau) / sqrt_lq;

    if (abs(costk) <= 1) {
        sintk = TMath::Sqrt(1. - costk * costk);
    } else {
        std::cout << "     rv2ln: costk > 1 " << costk << std::endl;
        sintk = 0;
    }

    mu = (E_h - p_l * costk - p_h * sintk * TMath::Cos(phi_rad - phi)) / M;
    factor = 1. + tau - mu;

    tails(tau, tm, mu);

    Double_t sfm[4];
    Double_t res = 0;

    if (ita == 1) {
        strf(0, 0, 0, sfm);
//        Double_t rmin = TMath::Power(10, -8);
//        Double_t rmax = (px2 - kMassC2) / factor;
//        simpux(rmin, rmax, 100, epsrr, podinl, res);
    } else if (ita == 2) {
        res = podinl((px2 - M * M) / factor) / factor;
    }

    return res * (Q2 / S_x + tau);
*/
    return 0;
}



void TRadCor::strf(Double_t tau, Double_t mu, Double_t R, Double_t (&sfm)[4])
{
    /*
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

    for (Int_t i = 0 ; i < 4 ; ++i) {
        sfm[i] = 0;
    }

    Double_t _Q2  = Q2 + R * tau;
    Double_t _t   = t - R * (tau - mu);
    Double_t _lq  = (S_x - R) * (S_x - R) + 4. * M * M * _Q2;
    Double_t _px2 = px2 - R * (1. + tau - mu);
    Double_t phq  = (m_h * m_h - _Q2 - _t) / 2.;

    Double_t _nu  = (S_x - R) / (2 * M);
    Double_t aks  = _Q2 / (2 * M) / _nu;
    Double_t zh   = E_h / _nu;
    Double_t _sq  = TMath::Sqrt(_Q2 + TMath::Power(_nu, 2));
    Double_t _p_l = (E_h * _nu - phq) / _sq;
    Double_t _pt2 = TMath::Power(p_h, 2) - TMath::Power(_p_l, 2);

    const Double_t& e = kEpsMachine;

    Double_t eps_nu = TMath::Sqrt(
                TMath::Power(S_x * e / (2 * M), 2) +
                TMath::Power(R * e / (2 * M), 2) +
                TMath::Power((S_x - R) / (2 * M) * e, 2));

    Double_t epst = TMath::Sqrt(
                TMath::Power((t * e), 2) +
                TMath::Power(((tau - mu) * R * e), 2) +
                TMath::Power((R * tau * e), 2) +
                TMath::Power((R * mu * e), 2));

    Double_t epsphq = TMath::Sqrt(
                TMath::Power((2 * m_h * m_h * e), 2) +
                TMath::Power((_Q2 * e), 2) +
                TMath::Power(epst, 2)
                ) / 2;

    Double_t epsq = 1. / 2. / TMath::Sqrt(_Q2 + _nu * _nu) *
                TMath::Sqrt(TMath::Power((_Q2 * e), 2) +
                        TMath::Power(2 * _nu * _nu * eps_nu, 2));

    Double_t epspl = TMath::Sqrt(
                TMath::Power(E_h * e * _nu / _sq, 2) +
                TMath::Power(E_h * eps_nu / _sq, 2) +
                TMath::Power(epsphq / _sq, 2) +
                    TMath::Power((E_h * _nu - phq) /
                        TMath::Power(_sq, 2) * epsq, 2));

    Double_t epspt2 = 2. * TMath::Sqrt(
                TMath::Power(p_h, 4) * TMath::Power(e, 2) +
                TMath::Power(_p_l, 4) * TMath::Power(epspl, 2));

    if (_pt2 < 0 && (_pt2 * _pt2 - epspt2 * epspt2) > 0) return;


    Double_t a = S / (2 * M) * (S / (2 * M) - nu) * _Q2 / Q2;
    Double_t tlde = (_nu + TMath::Sqrt(TMath::Power(_nu, 2) + 4 * a)) / 2;
    Double_t tldy = _nu / tlde;

    Double_t H1z, H2z, H3z, H4z;

    HapradUtils::SemiInclusiveModel(_Q2, aks, tldy, zh, _pt2, _px2,
                                    _p_l, H1z, H2z, H3z, H4z);

    // Including kinematic shift due to photon emission
    Double_t _sqrt_lq = TMath::Sqrt(_lq);

    Double_t aa = (S_x - R) * (zh - 2 * M * _p_l / _sqrt_lq) / 2. / (M * M);

    Double_t h1 = zh * (2. * _Q2 * aks * H1z - _pt2 * H4z) /
                M / _Q2 / aks;

    Double_t h2 = 2 * zh * (2 * aks * (aks * H2z + aa * H3z) * _lq +
            H4z * (aa * aa * _lq - 2 * _pt2 * _Q2)) /
                    M / aks / _Q2 / _lq;

    Double_t h3 = 2 * H4z * zh / M / _Q2 / aks;

    Double_t h4 = -2 * zh * (aks * H3z + aa * H4z) / M / _Q2 / aks;

    sfm[0] = h1;
    sfm[1] = h2;
    sfm[2] = h3;
    sfm[3] = h4;
*/
}
