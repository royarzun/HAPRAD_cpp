#ifndef TRADCOR_H
#define TRADCOR_H

#include "TROOT.h"


struct Epsilon {
    double epsphir, epstau, epsrr;

    Epsilon(): epsphir(0.01), epstau(0.001), epsrr(0.001) {}
};


struct Sxy {
    double s, x, sx, sxp, y, ym, w2;
    double als, alx, alm, aly, anu;
    double sqls, sqlx, sqly, sqlm, allm, an;
    double tamin, tamax, xs, ys;

    Sxy() :
        s(0.), x(0.), sx(0.), sxp(0.), y(0.), ym(0.), w2(0.),
        als(0.), alx(0.), alm(0.), aly(0.), anu(0.),
        sqls(0.), sqlx(0.), sqly(0.), sqlm(0.), allm(0.), an(0.),
        tamin(0.), tamax(0.), xs(0.), ys(0.) {}
};


struct Phi {
    double phirad, zdif, tdif, phidif;
    double p22, ehad, pph, plh, pth;
    double vv10, vv20, phk12, phkp, phkm, tdmin;
    int ilep;

    Phi():
        phirad(0.), zdif(0.), tdif(0.), phidif(0.),
        p22(0.), ehad(0.), pph(0.), plh(0.), pth(0.),
        vv10(0.), vv20(0.), phk12(0.), phkp(0.), phkm(0.), tdmin(0.),
        ilep(1) {}
};


struct Tail {
    double un, pl;
    int ita, isf1, isf2, isf3, ipol, iphi_rad, iphi_had;

    Tail():
        un(0.), pl(0.),
        ita(0), isf1(0), isf2(0), isf3(0),
        ipol(0), iphi_rad(0), iphi_had(0) {}
};



class TRadCor {
public:
    TRadCor();
    TRadCor(Double_t Ebeam, Double_t x, Double_t Q2, Double_t z,
            Double_t pt, Double_t phi, Double_t maxMx2);
    ~TRadCor();

    Double_t GetRCFactor(void);
    Double_t GetRCFactor(Double_t Ebeam, Double_t x, Double_t Q2, Double_t z,
                         Double_t pt, Double_t phi, Double_t maxMx2);

    void SetParameters(Double_t Ebeam, Double_t x, Double_t Q2, Double_t z,
                       Double_t pt, Double_t phi, Double_t maxMx2);
    void SetX(Double_t x);
    void SetQ2(Double_t Q2);
    void SetZ(Double_t z);
    void SetPt(Double_t pt);
    void SetPhi(Double_t phi);
    void SetEbeam(Double_t Ebeam);
    void SetMaxMx(Double_t maxMx2);

    void RegisteredLepton(Int_t type = 1);
    void IntegratePhiRad(Int_t type = 0);
    void IntegratePhiHad(Int_t type = 0);
    void SetPolarization(Int_t type = 0);

private:
    void Setup(void);
    void Haprad(void);
    void Conkin(void);
    void SPhiH(void);
    void Deltas(const Double_t massLepton2);
    Double_t VacPol(void);
    void Bornin(void);
    void BorninTest(Double_t& sib);
    void qqt(Double_t& tai);

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
    Double_t fDeltaInf;
    Double_t fDelta;
    Double_t fTail;
    Double_t fTai[3];

    Epsilon _eps;
    Phi _phi;
    Sxy _Sxy;
    Tail _tail;
};

#endif
