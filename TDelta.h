#include "TROOT.h"

class TGlobalConfig;
class TLorentzInvariants;
class THadronKinematics;


class TDelta {
public:
    TDelta();
    ~TDelta();

    Double_t    VR(void)  const { return fVR; };
    Double_t    Inf(void) const { return fInf; };
    Double_t    Vac(void) const { return fVac; };

    void        Evaluate(const TLorentzInvariants& inv,
                         const THadronKinematics& hkin);

    void        GlobalConfig(const TGlobalConfig* config) { fConfig = config; };

private:
    Double_t    VacPol(const Double_t Q2);

    const TGlobalConfig* fConfig;

    Double_t    fVR;
    Double_t    fInf;
    Double_t    fVac;
};
