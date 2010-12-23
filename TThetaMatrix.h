#ifndef TTHETAMATRIX_H
#define TTHETAMATRIX_H

#include "TMatrixD.h"

class TRadCor;
class TGlobalConfig;
class TKinematicalVariables;
class TLorentzInvariants;
class THadronKinematics;


class TThetaMatrix : public TMatrixD {
public:
    TThetaMatrix(const TRadCor* rc);
    TThetaMatrix(Int_t rows, Int_t cols, const TRadCor* rc);
    ~TThetaMatrix();

    void        Evaluate(Double_t tau, Double_t mu,
                         Int_t ita, Double_t phi_k);

private:
    const TGlobalConfig*            fConfig;
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;
    const THadronKinematics*        fHadKin;
};

#endif
