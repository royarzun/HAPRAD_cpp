#ifndef TRADCOR_H
#define TRADCOR_H

#include "TROOT.h"


class TRadCor {
public:
    TRadCor();
    TRadCor(Double_t Ebeam, Double_t x, Double_t q2, Double_t z,
            Double_t pt, Double_t phi, Double_t maxMx2);
    ~TRadCor();

    Double_t GetRCFactor(void);
    Double_t GetRCFactor(Double_t Ebeam, Double_t x, Double_t q2, Double_t z,
                         Double_t pt, Double_t phi, Double_t maxMx2);

    void SetParameters(Double_t Ebeam, Double_t x, Double_t q2, Double_t z,
                       Double_t pt, Double_t phi, Double_t maxMx2);
    void SetX(Double_t x);
    void SetQ2(Double_t q2);
    void SetZ(Double_t z);
    void SetPt(Double_t pt);
    void SetPhi(Double_t phi);
    void SetEbeam(Double_t Ebeam);
    void SetMaxMx(Double_t maxMx2);

private:
    void Setup(void);
    void Haprad(void);

    Double_t fRCFac;
    Double_t fEbeam;        // The beam energy
    Double_t fX;
    Double_t fY;
    Double_t fZ;
    Double_t fPt;
    Double_t fPhi;

    Double_t fMaxMx2;       // The maximum allowable amount of missing mass squared
    Double_t fMx2;          // Missing mass squared

    Double_t fSib;
    Double_t fSig;
    Double_t fDelta;
    Double_t fTail;
};

#endif
