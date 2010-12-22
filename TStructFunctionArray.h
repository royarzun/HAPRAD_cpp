#ifndef TSTRUCTFUNCTIONARRAY_H
#define TSTRUCTFUNCTIONARRAY_H

#include "TArrayD.h"

class TRadCor;
class TKinematicalVariables;
class TLorentzInvariants;
class THadronKinematics;


class TStructFunctionArray : public TArrayD {
public:
    TStructFunctionArray(const TRadCor* rc);
    TStructFunctionArray(const TRadCor* rc, Int_t n);
    ~TStructFunctionArray();

    void        Evaluate(Double_t tau, Double_t mu, Double_t R);

    Double_t    GetH1() const { return fArray[0]; };
    Double_t    GetH2() const { return fArray[1]; };
    Double_t    GetH3() const { return fArray[2]; };
    Double_t    GetH4() const { return fArray[3]; };
    Double_t    GetH(Int_t i) const { return GetAt(i); };

private:
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;
    const THadronKinematics*        fHadKin;
};

#endif
