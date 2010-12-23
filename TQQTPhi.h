#ifndef TQQTPHI_H
#define TQQTPHI_H

#include "Math/IFunction.h"

class TRadCor;
class TKinematicalVariables;
class TLorentzInvariants;


class TQQTPhi : public ROOT::Math::IBaseFunctionOneDim {
public:
    TQQTPhi(const TRadCor* rc);
    ~TQQTPhi();

    virtual IBaseFunctionOneDim* Clone() const;

private:
    virtual double DoEval(double phi) const;

    const TRadCor*                  fRC;
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;

    double fTauMax;
    double fTauMin;

    double fTauArray[6];
};

#endif
