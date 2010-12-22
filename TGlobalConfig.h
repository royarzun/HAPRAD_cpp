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

    void        SetPolarization(Int_t type) { fPolType = type; };
    void        SetLepton(Int_t type) { fLepton = type; };
    void        SetIntegrationPhiRad(Int_t type) { fPhiHad = type; };
    void        SetIntegrationPhiHad(Int_t type) { fPhiRad = type; };

private:
    Int_t       fPolType;
    Int_t       fLepton;
    Int_t       fPhiHad;
    Int_t       fPhiRad;
};

#endif
