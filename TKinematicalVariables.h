#include "TROOT.h"

class TKinematicalVariables {
public:
    TKinematicalVariables();
    TKinematicalVariables(Double_t x, Double_t Q2, Double_t z,
                          Double_t p_t, Double_t phi);
    ~TKinematicalVariables();

    Double_t    X(void)    const { return fX; };
    Double_t    Y(void)    const { return fY; };
    Double_t    Z(void)    const { return fZ; };
    Double_t    T(void)    const { return fT; };
    Double_t    PhiH(void) const { return fPhiH; };

    void        SetX(Double_t x) { fX = x; };
    void        SetY(Double_t y) { fY = y; };
    void        SetZ(Double_t z) { fZ = z; };
    void        SetT(Double_t t) { fT = t; };
    void        SetPhiH(Double_t phi) { fPhiH = phi; };

    void        SetAll(Double_t x, Double_t Q2, Double_t z,
                       Double_t p_t, Double_t phi);


private:
    Double_t    fX;
    Double_t    fY;
    Double_t    fZ;
    Double_t    fT;
    Double_t    fPhiH;
};
