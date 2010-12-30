#ifndef TSFFUN_H
#define TSFFUN_H

#include "TArrayD.h"

class TRadCor;
class TKinematicalVariables;
class TLorentzInvariants;
class THadronKinematics;

class TSffun : public TArrayD {
public:
    TSffun(const TRadCor *rc);
    ~TSffun();

    void        Evaluate(Double_t Q2, Double_t w2, Double_t t);

    Double_t    GetSfm0() const {
        return fArray[0];
    };
    Double_t    GetSfm1() const {
        return fArray[1];
    };
    Double_t    GetSfm2() const {
        return fArray[2];
    };
    Double_t    GetSfm3() const {
        return fArray[3];
    };
    Double_t    GetH(Int_t i) const {
        return GetAt(i);
    };
private:
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;
    const THadronKinematics*        fHadKin;
};


#endif
