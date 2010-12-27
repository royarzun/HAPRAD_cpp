#ifndef TRV2LN_H
#define TRV2LN_H

#include "Math/IFunction.h"
#include "TStructFunctionArray.h"

class TRadCor;
class TGlobalConfig;
class TKinematicalVariables;
class TLorentzInvariants;
class THadronKinematics;


class TRV2LN : public ROOT::Math::IBaseFunctionOneDim {
public:
    TRV2LN(const TRadCor* rc, double phi_k);
    ~TRV2LN();

    virtual IBaseFunctionOneDim* Clone() const;

private:
    virtual double DoEval(double tauln) const;

    const TRadCor*                  fRC;
    const TGlobalConfig*            fConfig;
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;
    const THadronKinematics*        fHadKin;

    TStructFunctionArray fH;

    double fPhiK;
    double fTauMax;
    double fTauMin;
};

#endif
