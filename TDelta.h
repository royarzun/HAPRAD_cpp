#ifndef TDELTA_H
#define TDELTA_H

#include "TROOT.h"

class TRadCor;
class TGlobalConfig;
class TKinematicalVariables;
class TLorentzInvariants;
class THadronKinematics;


class TDelta {
public:
    TDelta(const TRadCor* rc);
    ~TDelta();

    Double_t    VR(void)  const { return fVR; };
    Double_t    Inf(void) const { return fInf; };
    Double_t    Vac(void) const { return fVac; };

    void        Evaluate(void);

private:
    Double_t    VacPol(const Double_t Q2);

    const TGlobalConfig*            fConfig;
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;
    const THadronKinematics*        fHadKin;

    Double_t    fVR;
    Double_t    fInf;
    Double_t    fVac;
};

#endif
