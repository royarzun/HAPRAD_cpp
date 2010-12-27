#ifndef TGLOBALCONFIG_H
#define TGLOBALCONFIG_H

#include "TROOT.h"


class TGlobalConfig {
public:
    TGlobalConfig();

    Int_t       PolarizationType(void) const { return fPolType; };
    Int_t       LeptonType(void) const { return fLepton; };
    Int_t       IntegratePhiHad(void) const { return fPhiHad; };
    Int_t       IntegratePhiRad(void) const { return fPhiRad; };

    Int_t       EpsPhiR(void) const { return fEpsPhiR; };
    Int_t       EpsTau(void) const { return fEpsTau; };
    Int_t       EpsRR(void) const { return fEpsRR; };


    void        SetPolarization(Int_t type) { fPolType = type; };
    void        SetLepton(Int_t type) { fLepton = type; };
    void        SetIntegrationPhiRad(Int_t type) { fPhiHad = type; };
    void        SetIntegrationPhiHad(Int_t type) { fPhiRad = type; };

    void       SetEpsPhiR(Double_t value) { fEpsPhiR = value; };
    void       SetEpsTau(Double_t value) { fEpsTau = value; };
    void       SetEpsRR(Double_t value) { fEpsRR = value; };

private:
    Int_t       fPolType;
    Int_t       fLepton;
    Int_t       fPhiHad;
    Int_t       fPhiRad;

    // Epsilons
    Double_t    fEpsPhiR;
    Double_t    fEpsTau;
    Double_t    fEpsRR;
};

#endif
