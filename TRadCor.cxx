#include "TRadCor.h"
#include "THapradUtils.h"
#include "haprad_constants.h"
#include <iostream>
#include <iomanip>


TRadCor::TRadCor()
: E(0), maxMx2(0.), Mx2(0.),
  x(0.), y_i(0.), z(0.), t_i(0.), phi(0.),
  sib(0.), sig(0.), delta(0.), tail(0.),
  M(kMassProton), m(kMassElectron), m_h(kMassDetectedHadron),
  eps_phir(0.01), eps_tau(0.001), eps_rr(0.001),
  polType(0), int_phi_rad(0), int_phi_had(0)
{
    // Default constructor
}



TRadCor::TRadCor(Double_t E, Double_t x, Double_t Q2, Double_t z,
                 Double_t p_t, Double_t phi, Double_t maxMx2)
: E(0), maxMx2(0.), Mx2(0.),
  x(0.), y_i(0.), z(0.), t_i(0.), phi(0.),
  sib(0.), sig(0.), delta(0.), tail(0.),
  M(kMassProton), m(kMassElectron), m_h(kMassDetectedHadron),
  eps_phir(0.01), eps_tau(0.001), eps_rr(0.001),
  polType(0), int_phi_rad(0), int_phi_had(0)
{
    // Normal constructor for a radiative correction object
    //
    // The E parameter is the energy of the beam, the x, Q2, z, p_t and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount
    // of missing mass.

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

    this->E   = E;

    this->x   = x;
    this->y_i = -Q2;
    this->z   = z;
    this->t_i = p_t;
    this->phi = phi / kRadianDeg;

    this->maxMx2 = maxMx2;

    Setup();
}



void TRadCor::SetEbeam(Double_t E)
{
    this->E = E;
}



void TRadCor::SetX(Double_t x)
{
    this->x = x;
    Setup();
}



void TRadCor::SetQ2(Double_t Q2)
{
    this->y_i = -Q2;
    Setup();
}



void TRadCor::SetZ(Double_t z)
{
    this->z = z;
    Setup();
}



void TRadCor::SetPt(Double_t p_t)
{
    this->t_i = p_t;
    Setup();
}



void TRadCor::SetPhi(Double_t phi)
{
    this->phi = phi / kRadianDeg;
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

    switch(type) {
        case 1:
            m = kMassElectron;
            break;
        case 2:
            m = kMassMuon;
            break;
    }
}



void TRadCor::IntegratePhiRad(Int_t type)
{
    // Set whether to integrate over phi_{rad} (1) or to approximate (0).
    //
    // Default is 0;

    int_phi_rad = type;
}



void TRadCor::IntegratePhiHad(Int_t type)
{
    // Set whether to integrate over phi_{had} (1) or not (0).
    //
    // Default is 0;

    int_phi_had = type;
}



void TRadCor::SetPolarization(Int_t type)
{
    // Set the type of the polarization:
    //
    //     1 -- long; 2 -- tran; 0 -- unpol.
    //
    // Default is 0;

    polType = type;
}



void TRadCor::Setup(void)
{
    // Calculate the missing mass squared

    Double_t S_x = - y_i / x;
    Mx2 = M * M + S_x * (1.0 - z) + t_i;
}



Double_t TRadCor::GetRCFactor(void)
{
    // Get the radiative correction factor. You must set the parameters before
    // using this method.

    if (Mx2 > maxMx2) {
        Haprad();
        rc = sig / sib;
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
    isf1 = 1;
    isf2 = 4;
    isf3 = 1;

    pl = 0.;

    S = 2. * M * E;

    if (y_i >= 0.) {
        y = y_i;
        Q2 = S * x * y;
    } else {
        Q2 = - y_i;
        y = Q2 / (S * x);
    }

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "S      " << std::setw(20) << std::setprecision(10) << S  << std::endl;
    std::cout << "y      " << std::setw(20) << std::setprecision(10) << y << std::endl;
    std::cout << "Q^2    " << std::setw(20) << std::setprecision(10) << Q2  << std::endl;
#endif

    Double_t y_max = 1. / (1. + M * M * x / S);
    Double_t y_min = (kMassC2 - M * M) / (S * (1. - x));

    if (y > y_max || y < y_min || x > 1. || x < 0.) {
        std::cout << " Warning! Wrong kinematics! Skip the point!"
                  << std::endl
                  << " y = " << y << std::endl
                  << " x = " << x << std::endl;
        return;
    }

    Conkin();

    E_h = nu * z;
#ifdef DEBUG
    std::cout << "Eh     " << std::setw(20) << std::setprecision(10) << E_h  << std::endl;
#endif

    Double_t sqnuq = TMath::Sqrt(nu * nu + Q2);

    if (E_h < m_h) {
        std::cout << " Warning! Wrong kinematics! Skip the point!"
                  << std::endl
                  << " E_h =" << E_h
                  << std::endl;
        return;
    }

    p_h = TMath::Sqrt(E_h * E_h - m_h * m_h);
#ifdef DEBUG
    std::cout << "Ph     " << std::setw(20) << std::setprecision(10) << p_h  << std::endl;
#endif

    t = t_i;

    if (t >= 0.) {
        p_t = t;

        if (p_h < p_t) {
            std::cout << " Warning! Wrong kinematics! Skip the point!"
                      << std::endl
                      << " p_h = " << p_h << std::endl
                      << " p_t = " << p_t << std::endl;
            return;
        }

        p_l = TMath::Sqrt(p_h * p_h - p_t * p_t);

        if (p_h > p_t) {
            Double_t sqrt_lq = TMath::Sqrt(TMath::Max(0., lambda_q));
            N = N * sqrt_lq / 2. / M / p_l;
        } else {
            N = 0.;
        }

        t = m_h * m_h - Q2 + 2. * (sqnuq * p_l - nu * E_h);

        std::cout << "p_l: " << p_l << "\t" << t << "\t"
                  << p_l - t + Q2 - m_h * m_h + 2. * nu * E_h / 2. / sqnuq
                  << std::endl;
    } else {
        p_l = (t + Q2 - m_h * m_h + 2. * nu * E_h) / 2. / sqnuq;

        if (p_h < TMath::Abs(p_l)) {
            Double_t eps1, eps2, eps3, eps4, eps5, sum;

            eps1 = t * kEpsMachine / sqnuq;
            eps2 = 2. * m_h * m_h  * kEpsMachine / sqnuq;
            eps3 = 2. * nu * E_h * kEpsMachine / sqnuq;
            eps4 = t + Q2 - m_h * m_h + 2. * nu * E_h;
            eps5 = eps4 / sqnuq * kEpsMachine;

            sum = eps1 * eps1 + eps2 * eps2 + 2. * eps3 * eps3 + eps5 * eps5;

            Double_t epspl  = TMath::Sqrt(sum) / 2.;
            Double_t calEps = p_h - TMath::Abs(p_l);
            if (TMath::Abs(calEps) > epspl) {
               std::cout << " Warning! Wrong kinematics! Skeep the point!"
                         << std::endl << " p_h  = " << p_h
                         << std::endl << " p_l  = " << p_l
                         << "\t"      << calEps << std::endl;
               return;
            } else {
               std::cout << "Zero p_t! " << p_l
                         << "\t"         << calEps
                         << "\t"         << epspl << std::endl;
               p_l = TMath::Sign(1., p_l) * p_h;
            }
        }
        p_t = TMath::Sqrt(p_h * p_h - p_l * p_l);
    }

#ifdef DEBUG
    std::cout << "Pt     " << std::setw(20) << std::setprecision(10) << p_t  << std::endl;
    std::cout << "Pl     " << std::setw(20) << std::setprecision(10) << p_l  << std::endl;
    std::cout << "tdif   " << std::setw(20) << std::setprecision(10) << t << std::endl;
#endif

    Double_t px2_max = W2 - m_h * m_h;
    px2 = M * M + S_x * (1. - z) + t;

    if (px2 < kMassC2 || px2 > px2_max) {
        std::cout << " Warning! Wrong kinematics! Skeep the point!" << std::endl
                  << " p_x^2     = " << px2 << std::endl
                  << " mc2       = " << kMassC2 << std::endl
                  << " p_x^2 max = " << px2_max << std::endl;
        return;
    }

    t_min = m_h * m_h - Q2 + 2. * (sqnuq * p_h - nu * E_h);

    Double_t tdmax;
    tdmax = m_h * m_h - Q2 + 2. * (- sqnuq * p_h - nu * E_h);

    if ((t - t_min) > kEpsMachine || t < tdmax) {
        std::cout << " Warning! Wrong kinematics! Skeep the point!" << std::endl
                  << " t     = " << t << std::endl
                  << " t_max = " << tdmax << std::endl
                  << " t_min = " << t_min << " - " << t - t_min << std::endl;
        return;
    }

    SPhiH();
}



void TRadCor::Conkin()
{
    X   = S * (1. - y);
    S_x = S - X;
    S_p = S + X;
    W2  = S - X - Q2 + M * M;

    lambda_s = S * S - 4. * m * m * M * M;
    lambda_x = X * X - 4. * m * m * M * M;
    lambda_m = Q2 * Q2 + 4. * m * m * Q2;
    lambda_q = TMath::Power(S_x, 2) + 4. * M * M * Q2;

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "X      " << std::setw(20) << std::setprecision(10) << X        << std::endl;
    std::cout << "S_x    " << std::setw(20) << std::setprecision(10) << S_x      << std::endl;
    std::cout << "S_p    " << std::setw(20) << std::setprecision(10) << S_p      << std::endl;
    std::cout << "W2     " << std::setw(20) << std::setprecision(10) << W2       << std::endl;
    std::cout << "l_s    " << std::setw(20) << std::setprecision(10) << lambda_s << std::endl;
    std::cout << "l_x    " << std::setw(20) << std::setprecision(10) << lambda_x << std::endl;
    std::cout << "l_m    " << std::setw(20) << std::setprecision(10) << lambda_m << std::endl;
    std::cout << "l_q    " << std::setw(20) << std::setprecision(10) << lambda_q << std::endl;
#endif

    if (lambda_s < 0.) std::cout << " Conkin: lambda_s < 0 " << std::endl;
    if (lambda_x < 0.) std::cout << " Conkin: lambda_x < 0 " << std::endl;
    if (lambda_q < 0.) std::cout << " Conkin: lambda_q < 0 " << std::endl;
    if (lambda_m < 0.) std::cout << " Conkin: lambda_m < 0 " << std::endl;


    Double_t sqrt_lq, sqrt_lm;
    sqrt_lq = TMath::Sqrt(TMath::Max(0., lambda_q));
    sqrt_lm = TMath::Sqrt(TMath::Max(0., lambda_m));

    nu  = S_x / (2. * M);
    N   = kPi * kAlpha * kAlpha * y * S_x * M / 2. / sqrt_lq * kBarn;

#ifdef DEBUG
    std::cout << "nu     " << std::setw(20) << std::setprecision(10) << nu << std::endl;
    std::cout << "N      " << std::setw(20) << std::setprecision(10) << N << std::endl;
#endif

    tau_max = (S_x + sqrt_lq) / (2. * M * M);
    tau_min = - Q2 / (M * M) / tau_max;
}



void TRadCor::SPhiH(void)
{
    Double_t sqrt_ls, sqrt_lx, sqrt_lq;
    Double_t costs, costx, sints, sintx;
    Double_t lambda;

    sqrt_ls = TMath::Sqrt(TMath::Max(0., lambda_s));
    sqrt_lx = TMath::Sqrt(TMath::Max(0., lambda_x));
    sqrt_lq = TMath::Sqrt(TMath::Max(0., lambda_q));

    costs = (S * (S - X) + 2. * M * M * Q2) / sqrt_ls / sqrt_lq;
    costx = (X * (S - X) - 2. * M * M * Q2) / sqrt_lx / sqrt_lq;

    lambda = S * X * Q2 - M * M * Q2 * Q2 - m * m * lambda_q;

    if (lambda > 0) {
        sints = 2. * M * TMath::Sqrt(lambda) / sqrt_ls / sqrt_lq;
        sintx = 2. * M * TMath::Sqrt(lambda) / sqrt_lx / sqrt_lq;
    } else {
        std::cout << "sphi: sints = NaN " << lambda << std::endl;
        std::cout << "sphi: sintx = NaN " << lambda << std::endl;
        sints = 0.;
        sintx = 0.;
    }

    Double_t v1, v2;
    v1 = costs * p_l + sints * p_t * TMath::Cos(phi);
    v2 = costx * p_l + sintx * p_t * TMath::Cos(phi);

    V_1 = (S * E_h - sqrt_ls * v1) / M;
    V_2 = (X * E_h - sqrt_lx * v2) / M;

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "V1     " << std::setw(20) << std::setprecision(10) << V_1  << std::endl;
    std::cout << "V2     " << std::setw(20) << std::setprecision(10) << V_2  << std::endl;
#endif

    Deltas();

    Double_t sibt;
    Int_t it_end = 3;

    if (polType == 0) {
        it_end = 1;
    }

    for (Int_t i = 1; i <= it_end; ++i) {
        for (ita = 1; ita <= 2; ++ita) {
            std::cout << "********** ita: " << ita
                      << " *********" << std::endl;
            if (ita == 1) {
                sib = Bornin();
                BorninTest(sibt);
                std::cout << "sib1" << sib << std::endl;
                std::cout << "sibt" << sibt << std::endl;
                if (sib == 0.0) {
                    tai[1] = 0.;
                    continue;
                }
            }
            qqt(tai[ita]);
            std::cout << "tai[" << ita
                      << "]\t"  << tai[ita] << std::endl;
        }
         del_inf = 0.;
         Double_t extai1 = TMath::Exp(kAlpha / kPi * del_inf);
         sig = sib * extai1 * (1. + kAlpha / kPi * (delta - del_inf)) +
                    tai[1] + tai[2];
    }
}



void TRadCor::Deltas(void)
{
    Double_t delta_vac = VacPol();
    Double_t S_ = S - Q2 - V_1;
    Double_t X_ = X + Q2 - V_2;

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "d_vac  " << std::setw(20) << std::setprecision(10) << delta_vac << std::endl;
    std::cout << "S'     " << std::setw(20) << std::setprecision(10) << S_  << std::endl;
    std::cout << "X'     " << std::setw(20) << std::setprecision(10) << X_  << std::endl;
#endif

    Double_t lambda_s_ = TMath::Power(S_, 2) - 4. * px2 * m * m;
    Double_t lambda_x_ = TMath::Power(X_, 2) - 4. * px2 * m * m;
    if (lambda_s_ < 0.)
        std::cout << "deltas: lambda_s' < 0 " << lambda_s_ << std::endl;
    if (lambda_x_ < 0.)
        std::cout << "deltas: lambda_x' < 0 " << lambda_x_ << std::endl;

    Double_t l_m = TMath::Log(Q2 / (m * m));
    Double_t Li_2 = HapradUtils::fspen(1. - px2 * Q2 / (S_ * X_));
    Double_t delta_VR = 1.5 * l_m - 2. -
            0.5 * TMath::Power(TMath::Log(X_ / S_),2) + Li_2 - kPi * kPi / 6.;

    del_inf = (l_m - 1.) * TMath::Log(TMath::Power(px2 - kMassC2, 2) / S_ / X_);
    delta   = del_inf + delta_vac + delta_VR;
}



Double_t TRadCor::VacPol(void)
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



Double_t TRadCor::Bornin(void)
{
    Double_t sfm[8];
    Double_t tm[isf2];

    strf(0., 0., 0, sfm);

    tm[0] = Q2;
    tm[1] = (S * X - M * M * Q2) / 2.;
    tm[2] = (V_1 * V_2 - m_h * m_h * Q2) / 2.;
    tm[3] = (S * V_2 + X * V_1 - z * S_x * Q2) / 2.;

#ifdef DEBUG
        std::cout << "    BORNIN      " << std::endl;
#endif

    Double_t ssum = 0.;
    for (Int_t i = isf1; i <= isf2; i += isf3) {
        ssum = ssum + tm[i] * sfm[i];
#ifdef DEBUG
        std::cout.setf(std::ios::fixed);
        std::cout << "    i      " << std::setw(20) << std::setprecision(10) << i      << std::endl;
        std::cout << "    tm[i]  " << std::setw(20) << std::setprecision(10) << tm[i]  << std::endl;
        std::cout << "    sfm[i] " << std::setw(20) << std::setprecision(10) << sfm[i] << std::endl;
        std::cout << "    ssum   " << std::setw(20) << std::setprecision(10) << ssum   << std::endl;
#endif
    }

    return ssum * N / Q2 / Q2 * 2.;
}



void TRadCor::BorninTest(Double_t& sib)
{

}



void TRadCor::qqt(Double_t& tai)
{

}



void TRadCor::strf(Double_t tau, Double_t d2kvir, Double_t R, Double_t (&sfm)[8])
{
}
