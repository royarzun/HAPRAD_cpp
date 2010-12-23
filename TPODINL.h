#ifndef TPODINL_H
#define TPODINL_H

#include "Math/IFunction.h"

class TRadCor;
class TLorentzInvariants;
class TStructFunctionArray;
class TThetaMatrix;


class TPODINL : public ROOT::Math::IBaseFunctionOneDim {
public:
    TPODINL(const TRadCor* rc, double tau, double mu,
            const TStructFunctionArray& H0, const TThetaMatrix& theha);
    ~TPODINL();

    virtual IBaseFunctionOneDim* Clone() const;

private:
    virtual double DoEval(double R) const;

    const TRadCor*                  fRC;
    const TLorentzInvariants*       fInv;

    double fTau;
    double fMu;
    const TStructFunctionArray& fH0;
    const TThetaMatrix& fTheta;
};

#endif
