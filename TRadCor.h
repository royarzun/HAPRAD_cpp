#ifndef TRADCOR_H
#define TRADCOR_H

#include "TROOT.h"


class TRadCor {
public:
    TRadCor();
    TRadCor(Double_t E, Double_t x, Double_t Q2, Double_t z,
            Double_t p_t, Double_t phi, Double_t maxMx2);
    ~TRadCor();

    Double_t    GetRCFactor(void);
    Double_t    GetRCFactor(Double_t E, Double_t x, Double_t Q2, Double_t z,
                            Double_t p_t, Double_t phi, Double_t maxMx2);

    void        SetParameters(Double_t E, Double_t x, Double_t Q2, Double_t z,
                              Double_t p_t, Double_t phi, Double_t maxMx2);
    void        SetX(Double_t x);
    void        SetQ2(Double_t Q2);
    void        SetZ(Double_t z);
    void        SetPt(Double_t p_t);
    void        SetPhi(Double_t phi);
    void        SetEbeam(Double_t E);
    void        SetMaxMx(Double_t maxMx2);

    void        RegisteredLepton(Int_t type = 1);
    void        IntegratePhiRad(Int_t type = 0);
    void        IntegratePhiHad(Int_t type = 0);
    void        SetPolarization(Int_t type = 0);

private:
    void        Setup(void);
    void        Haprad(void);
    void        Conkin(void);
    void        SPhiH(void);
    void        Deltas(void);
    Double_t    VacPol(void);
    Double_t    Bornin(void);
    void        BorninTest(Double_t& sigma_born);
    void        qqt(Double_t& tai);
    Double_t    qqtphi(Double_t phi);
    Double_t    rv2ln(Double_t tauln);
    Double_t    podinl(Double_t R);
    Double_t    tails(Double_t tau, Double_t (&theta)[4][3], Double_t mu);
    void        strf(Double_t tau, Double_t mu, Double_t R, Double_t (&sfm)[4]);


    Double_t     E;        // The beam energy
    Double_t     maxMx2;       // The maximum allowable amount of missing mass squared
    Double_t     Mx2;          // Missing mass squared

    //  Kinematic variables
    Double_t     x;
    Double_t     y;
    Double_t     y_i;
    Double_t     z;
    Double_t     t;
    Double_t     t_i;
    Double_t     t_min;
    Double_t     phi;

    // Results
    Double_t     rc;
    Double_t     sigma_born;    // sigma_0
    Double_t     sig_obs;       // sigma_{obs}
    Double_t     del_inf;
    Double_t     delta;
    Double_t     tail;
    Double_t     tai[3];

    // Masses
    Double_t     M;
    Double_t     m;
    Double_t     m_h;

    // Lorentz invariants
    Double_t     S;
    Double_t     X;
    Double_t     Q2;
    Double_t     W2;
    Double_t     S_x;
    Double_t     S_p;
    Double_t     lambda_q;
    Double_t     lambda_s;
    Double_t     lambda_x;
    Double_t     V_1;
    Double_t     V_2;

    Double_t     lambda_m;

    // Detected hadron kinematics
    Double_t     E_h;
    Double_t     p_h;
    Double_t     p_t;
    Double_t     p_l;
    Double_t     nu;

    // Deltas
    Double_t     px2;

    // Integration
    Double_t     N;                // Normalization factor
    Double_t     pl;
    Int_t        ita;
    Double_t     tau_max;
    Double_t     tau_min;
    Double_t     phi_rad;

    // Epsilons
    Double_t     eps_phir;
    Double_t     eps_tau;
    Double_t     eps_rr;

    // Configurations
    Int_t        polType;
    Int_t        int_phi_rad;
    Int_t        int_phi_had;
};

#endif
