#ifndef TRADCOR_H
#define TRADCOR_H

#include "TGlobalConfig.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
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

    const TGlobalConfig*            GetConfig(void) const { return &fConfig; };
    const TKinematicalVariables*    GetKinematicalVariables(void) const { return &fKin; };
    const TLorentzInvariants*       GetLorentzInvariants(void) const { return &fInv; };
    const THadronKinematics*        GetHadronKinematics(void) const { return &fHadKin; };

private:
    void        Setup(void);
    void        Haprad(void);
    void        SPhiH(void);
    void        BorninTest(Double_t& sigma_born);
    void        qqt(Double_t& tai);

    TGlobalConfig           fConfig;
    TKinematicalVariables   fKin;
    TLorentzInvariants      fInv;
    THadronKinematics       fHadKin;

    Double_t     E;        // The beam energy
    Double_t     maxMx2;       // The maximum allowable amount of missing mass squared
    Double_t     Mx2;          // Missing mass squared

    //  Kinematic variables
    Double_t     t_min;

    // Results
    Double_t     rc;
    Double_t     sigma_born;    // sigma_0
    Double_t     sig_obs;       // sigma_{obs}
    Double_t     tail;
    Double_t     tai[3];

    // Masses
    Double_t     M;
    Double_t     m;
    Double_t     m_h;

    // Integration
    Double_t     N;                // Normalization factor
    Double_t     pl;
    Int_t        ita;
    Double_t     phi_rad;
};

#endif
