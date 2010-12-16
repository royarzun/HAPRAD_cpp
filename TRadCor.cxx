#include "TRadCor.h"
#include "THapradUtils.h"
#include "haprad_constants.h"
#include <iostream>
#include <iomanip>


TRadCor::TRadCor()
: fEbeam(0),
  fX(0), fY(0), fZ(0), fPt(0), fPhi(0),
  fMaxMx2(0), fMx2(0), fSib(0), fSig(0), fDelta(0), fTail(0),
  _eps(), _phi(), _Sxy(), _tail()
{
    // Default constructor
}



TRadCor::TRadCor(Double_t Ebeam, Double_t x, Double_t Q2, Double_t z,
                 Double_t pt, Double_t phi, Double_t maxMx2)
: fEbeam(Ebeam),
  fX(x), fY(-Q2), fZ(z), fPt(pt), fPhi(phi/kRadianDeg),
  fMaxMx2(0), fMx2(0), fSib(0), fSig(0), fDelta(0), fTail(0),
  _eps(), _phi(), _Sxy(), _tail()
{
    // Normal constructor for a radiative correction object
    //
    // The Ebeam parameter is the energy of the beam, the x, Q2, z, pt and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    Setup();
}



TRadCor::~TRadCor()
{
    // Default destructor
}



void TRadCor::SetParameters(Double_t Ebeam, Double_t x, Double_t Q2, Double_t z,
                            Double_t pt, Double_t phi, Double_t maxMx2)
{
    // Set the values for the variables used to calculate the radiative
    // correction.
    //
    // The Ebeam parameter is the energy of the beam, the x, Q2, z, pt and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    fEbeam = Ebeam;
    fX = x;
    fY = -Q2;
    fZ = z;
    fPt = pt;
    fPhi = phi / kRadianDeg;
    fMaxMx2 = maxMx2;

    Setup();
}



void TRadCor::SetEbeam(Double_t Ebeam)
{
    fEbeam = Ebeam;
}



void TRadCor::SetX(Double_t x)
{
    fX = x;
    Setup();
}



void TRadCor::SetQ2(Double_t Q2)
{
    fY = -Q2;
    Setup();
}



void TRadCor::SetZ(Double_t z)
{
    fZ = z;
    Setup();
}



void TRadCor::SetPt(Double_t pt)
{
    fPt = pt;
    Setup();
}



void TRadCor::SetPhi(Double_t phi)
{
    fPhi = phi / kRadianDeg;
}



void TRadCor::SetMaxMx(Double_t maxMx2)
{
    fMaxMx2 = maxMx2;
}



void TRadCor::RegisteredLepton(Int_t type)
{
    // Set the registere lepton: 1 -- electron; 2 -- muon.
    //
    // Default is 1.

    _phi.ilep = type;
}



void TRadCor::IntegratePhiRad(Int_t type)
{
    // Set whether to integrate over phi_{rad} (1) or to approximate (0).
    //
    // Default is 0;

    _tail.iphi_rad = type;
}



void TRadCor::IntegratePhiHad(Int_t type)
{
    // Set whether to integrate over phi_{had} (1) or not (0).
    //
    // Default is 0;

    _tail.iphi_had = type;
}



void TRadCor::SetPolarization(Int_t type)
{
    // Set the type of the polarization:
    //
    //     1 -- long; 2 -- tran; 0 -- unpol.
    //
    // Default is 0;

    _tail.ipol = type;
}



void TRadCor::Setup(void)
{
    // Calculate the missing mass squared

    Double_t Sx;

    Sx = - fY / fX;
    fMx2 = TMath::Power(kMassProton,2) + Sx * (1.0 - fZ) + fPt;
}



Double_t TRadCor::GetRCFactor(void)
{
    // Get the radiative correction factor. You must set the parameters before
    // using this method.

    if (fMx2 > fMaxMx2) {
        Haprad();
        fRCFac = fSig / fSib;
    } else {
        fRCFac = 0;
    }
    return fRCFac;
}



Double_t TRadCor::GetRCFactor(Double_t Ebeam, Double_t x, Double_t Q2, Double_t z,
                              Double_t pt, Double_t phi, Double_t maxMx2)
{
    // Get the radiative correction factor for the given parameters.
    //
    // The Ebeam parameter is the energy of the beam, the x, Q2, z, pt and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    SetParameters(Ebeam,x,Q2,z,pt,phi,maxMx2);
    return GetRCFactor();
}



void TRadCor::Haprad(void)
{
    _tail.isf1 = 1;
    _tail.isf2 = 4;
    _tail.isf3 = 1;

    _tail.un = 1.;
    _tail.pl = 0.;

    _Sxy.xs = fX;
    _phi.zdif = fZ;
    _phi.tdif = fPt;

    _Sxy.s = 2. * kMassProton * fEbeam;

    if (fY >= 0.) {
        _Sxy.ys = fY;
        _Sxy.y = _Sxy.s * _Sxy.xs * _Sxy.ys;
    } else {
        _Sxy.y = - fY;
        _Sxy.ys = _Sxy.y / (_Sxy.s * _Sxy.xs);
    }

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "S      " << std::setw(20) << std::setprecision(10) << _Sxy.s  << std::endl;
    std::cout << "y      " << std::setw(20) << std::setprecision(10) << _Sxy.ys << std::endl;
    std::cout << "Q^2    " << std::setw(20) << std::setprecision(10) << _Sxy.y  << std::endl;
#endif

    Double_t mp2 = TMath::Power(kMassProton, 2);
    Double_t yma = 1. / (1. + mp2 * _Sxy.xs / _Sxy.s);
    Double_t ymi = (kMassC2 - mp2) / (_Sxy.s * (1. - _Sxy.xs));

    if (_Sxy.ys > yma || _Sxy.ys < ymi || _Sxy.xs > 1. || _Sxy.xs < 0.) {
        std::cout << " Warning! Wrong kinematics!!!! skip the point!"
                  << std::endl
                  << " ys= " << _Sxy.ys << std::endl
                  << " xs= " << _Sxy.xs << std::endl;
        return;
    }

    Conkin();

    _phi.ehad = _Sxy.anu * _phi.zdif;
#ifdef DEBUG
    std::cout << "Eh     " << std::setw(20) << std::setprecision(10) << _phi.ehad  << std::endl;
#endif
    Double_t sqnuq = TMath::Sqrt(_Sxy.anu * _Sxy.anu + _Sxy.y);

    if (_phi.ehad < kMassDetectedHadron) {
        std::cout << " Warning! Wrong kinematics!!!! skeep the point!"
                  << std::endl
                  << " ehad =" << _phi.ehad
                  << std::endl;
        return;
    }

    Double_t mhh2 = TMath::Power(kMassDetectedHadron,2);

    _phi.pph = TMath::Sqrt(_phi.ehad * _phi.ehad - mhh2);
#ifdef DEBUG
    std::cout << "Ph     " << std::setw(20) << std::setprecision(10) << _phi.pph  << std::endl;
#endif

    if (_phi.tdif >= 0.) {
        _phi.pth = _phi.tdif;

        if (_phi.pph < _phi.pth) {
            std::cout << " Warning! Wrong kinematics!!!! skeep the point!"
                      << std::endl
                      << " pph =" << _phi.pph << std::endl
                      << " pth =" << _phi.pth << std::endl;
            return;
        }

        _phi.plh = TMath::Sqrt(_phi.pph * _phi.pph - _phi.pth * _phi.pth);

        if (_phi.pph > _phi.pth)
            _Sxy.an = _Sxy.an * _Sxy.sqly / 2. / kMassProton / _phi.plh;
        else
            _Sxy.an = 0.;

        _phi.tdif = mhh2 - _Sxy.y + 2. * (sqnuq * _phi.plh - _Sxy.anu * _phi.ehad);

        Double_t ptmp = _phi.tdif + _Sxy.y - mhh2 + 2. * _Sxy.anu * _phi.ehad;
        std::cout << "pl: " << _phi.plh
                  << "\t"   << _phi.tdif
                  << "\t"   << _phi.plh - ptmp / 2. / sqnuq
                  << std::endl;
    } else {
        Double_t ptmp = _phi.tdif + _Sxy.y - mhh2 + 2. * _Sxy.anu * _phi.ehad;
        _phi.plh = ptmp / 2. / sqnuq;

        if (_phi.pph < TMath::Abs(_phi.plh)) {
            Double_t eps1, eps2, eps3, eps4, eps5, sum;

            eps1 = _phi.tdif * kEpsMachine / sqnuq;
            eps2 = 2. * mhh2  * kEpsMachine / sqnuq;
            eps3 = 2. * _Sxy.anu * _phi.ehad * kEpsMachine / sqnuq;
            eps4 = _phi.tdif + _Sxy.y - mhh2 + 2. * _Sxy.anu * _phi.ehad;
            eps5 = eps4 / sqnuq * kEpsMachine;

            sum = eps1 * eps1 + eps2 * eps2 + 2. * eps3 * eps3 + eps5 * eps5;

            Double_t epspl  = TMath::Sqrt(sum) / 2.;
            Double_t calEps = _phi.pph - TMath::Abs(_phi.plh);
            if (TMath::Abs(calEps) > epspl) {
               std::cout << "Warning! Wrong kinematics! Skeep the point!"
                         << std::endl
                         << "pph  = " << _phi.pph
                         << std::endl
                         << "plh  = " << _phi.plh
                         << "\t"      << calEps
                         << std::endl;
               return;
            } else {
               std::cout << "Zero pt! " << _phi.plh
                         << "\t"        << calEps
                         << "\t"        << epspl
                         << std::endl;
               _phi.plh = TMath::Sign(1., _phi.plh) * _phi.pph;
            }
        }
        _phi.pth = TMath::Sqrt(_phi.pph * _phi.pph - _phi.plh * _phi.plh);
    }

#ifdef DEBUG
    std::cout << "Pt     " << std::setw(20) << std::setprecision(10) << _phi.pth  << std::endl;
    std::cout << "Pl     " << std::setw(20) << std::setprecision(10) << _phi.plh  << std::endl;
    std::cout << "tdif   " << std::setw(20) << std::setprecision(10) << _phi.tdif << std::endl;
#endif

    Double_t p22max = _Sxy.w2 - mhh2;
    _phi.p22 = mp2 + _Sxy.sx * (1. - _phi.zdif) + _phi.tdif;

    if (_phi.p22 < kMassC2 || _phi.p22 > p22max) {
        std::cout << "Warning! Wrong kinematics! Skeep the point!"
                  << std::endl
                  << "p22  = " << _phi.p22 << std::endl
                  << "amc2 = " << kMassC2 << std::endl
                  << "p22m = " << p22max << std::endl;
        return;
    }

    _phi.tdmin = mhh2 - _Sxy.y + 2. * (sqnuq * _phi.pph - _Sxy.anu * _phi.ehad);

    Double_t tdmax;
    tdmax = mhh2 - _Sxy.y + 2. * (- sqnuq * _phi.pph - _Sxy.anu * _phi.ehad);

    if ((_phi.tdif - _phi.tdmin) > kEpsMachine || _phi.tdif < tdmax) {
        std::cout << "Warning! Wrong kinematics! Skeep the point!"
                  << std::endl
                  << "tdif  = " << _phi.tdif << std::endl
                  << "tdmax = " << tdmax << std::endl
                  << "tdmin = " << _phi.tdmin
                  << " " << _phi.tdif - _phi.tdmin
                  << std::endl;
        return;
    }

    SPhiH();
}



void TRadCor::Conkin()
{
    Double_t massLepton;

    switch(_phi.ilep) {
        case 1:
            massLepton = kMassElectron;
            break;
        case 2:
            massLepton = kMassMuon;
            break;
        default:
            massLepton = kMassElectron;
    }

    Double_t ml2, mp2;
    ml2 = TMath::Power(massLepton,2);
    mp2 = TMath::Power(kMassProton,2);

    _Sxy.x   = _Sxy.s * (1 - _Sxy.ys);
    _Sxy.sx  = _Sxy.s - _Sxy.x;
    _Sxy.sxp = _Sxy.s + _Sxy.x;
    _Sxy.ym  = _Sxy.y + 2 * ml2;
    _Sxy.w2  = mp2 + _Sxy.s - _Sxy.y - _Sxy.x;
    _Sxy.als = _Sxy.s * _Sxy.s - 2 * ml2 * (2 * mp2);
    _Sxy.alx = _Sxy.x * _Sxy.x - 2 * ml2 * (2 * mp2);
    _Sxy.alm = _Sxy.y * _Sxy.y + 4 * ml2 * _Sxy.y;
    _Sxy.aly = TMath::Power(_Sxy.sx, 2) + 4 * mp2 * _Sxy.y;

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "x      " << std::setw(20) << std::setprecision(10) << _Sxy.x   << std::endl;
    std::cout << "sx     " << std::setw(20) << std::setprecision(10) << _Sxy.sx  << std::endl;
    std::cout << "sxp    " << std::setw(20) << std::setprecision(10) << _Sxy.sxp << std::endl;
    std::cout << "ym     " << std::setw(20) << std::setprecision(10) << _Sxy.ym  << std::endl;
    std::cout << "w2     " << std::setw(20) << std::setprecision(10) << _Sxy.w2  << std::endl;
    std::cout << "als    " << std::setw(20) << std::setprecision(10) << _Sxy.als << std::endl;
    std::cout << "alx    " << std::setw(20) << std::setprecision(10) << _Sxy.alx << std::endl;
    std::cout << "alm    " << std::setw(20) << std::setprecision(10) << _Sxy.alm << std::endl;
    std::cout << "aly    " << std::setw(20) << std::setprecision(10) << _Sxy.aly << std::endl;
#endif

    if (_Sxy.als < 0) std::cout << "Conkin: als < 0 " << std::endl;
    if (_Sxy.alx < 0) std::cout << "Conkin: alx < 0 " << std::endl;
    if (_Sxy.aly < 0) std::cout << "Conkin: aly < 0 " << std::endl;
    if (_Sxy.alm < 0) std::cout << "Conkin: alm < 0 " << std::endl;

    _Sxy.sqls = TMath::Sqrt(TMath::Max(0., _Sxy.als));
    _Sxy.sqlx = TMath::Sqrt(TMath::Max(0., _Sxy.alx));
    _Sxy.sqly = TMath::Sqrt(TMath::Max(0., _Sxy.aly));
    _Sxy.sqlm = TMath::Sqrt(TMath::Max(0., _Sxy.alm));

    Double_t lm = TMath::Log((_Sxy.sqlm + _Sxy.y) / (_Sxy.sqlm - _Sxy.y));
    _Sxy.allm = lm / _Sxy.sqlm;
    _Sxy.anu  = _Sxy.sx / (2 * kMassProton);
    _Sxy.an   = kPi * kAlpha * kAlpha * _Sxy.ys * _Sxy.sx *
                                    kMassProton / 2 / _Sxy.sqly * kBarn;

#ifdef DEBUG
    std::cout << "allm   " << std::setw(20) << std::setprecision(10) << _Sxy.allm << std::endl;
    std::cout << "anu    " << std::setw(20) << std::setprecision(10) << _Sxy.anu  << std::endl;
    std::cout << "an     " << std::setw(20) << std::setprecision(10) << _Sxy.an   << std::endl;
#endif

    _Sxy.tamax = (_Sxy.sx + _Sxy.sqly) / (2 * mp2);
    _Sxy.tamin = - _Sxy.y / mp2 / _Sxy.tamax;
}



void TRadCor::SPhiH(void)
{
    Double_t massLepton;

    switch(_phi.ilep) {
        case 1:
            massLepton = kMassElectron;
            break;
        case 2:
            massLepton = kMassMuon;
            break;
        default:
            massLepton = kMassElectron;
    }

    Double_t ml2, mp2;
    ml2 = TMath::Power(massLepton,2);
    mp2 = TMath::Power(kMassProton,2);

    Double_t costs, costx, sints, sintx;

    costs = (_Sxy.s * (_Sxy.s - _Sxy.x) + 2. * mp2 * _Sxy.y) / _Sxy.sqls / _Sxy.sqly;
    costx = (_Sxy.x * (_Sxy.s - _Sxy.x) - 2. * mp2 * _Sxy.y) / _Sxy.sqlx / _Sxy.sqly;

    Double_t lambda;
    lambda = _Sxy.s * _Sxy.x * _Sxy.y - mp2 * _Sxy.y * _Sxy.y - ml2 * _Sxy.aly;

    if (lambda > 0) {
        sints = 2. * kMassProton * TMath::Sqrt(lambda) / _Sxy.sqls / _Sxy.sqly;
        sintx = 2. * kMassProton * TMath::Sqrt(lambda) / _Sxy.sqlx / _Sxy.sqly;
    } else {
        std::cout << "sphi: sints = NaN " << lambda << std::endl;
        std::cout << "sphi: sintx = NaN " << lambda << std::endl;
        sints = 0.;
        sintx = 0.;
    }

    Double_t vv10, vv20;
    vv10 = costs * _phi.plh + sints * _phi.pth * TMath::Cos(fPhi);
    vv20 = costx * _phi.plh + sintx * _phi.pth * TMath::Cos(fPhi);

    _phi.vv10 = (_Sxy.s * _phi.ehad - _Sxy.sqls * vv10) / kMassProton;
    _phi.vv20 = (_Sxy.x * _phi.ehad - _Sxy.sqlx * vv20) / kMassProton;

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "V1     " << std::setw(20) << std::setprecision(10) << _phi.vv10  << std::endl;
    std::cout << "V2     " << std::setw(20) << std::setprecision(10) << _phi.vv20  << std::endl;
#endif

    Deltas(ml2);

    Double_t sibt;
    Double_t it_end = 3;

    if (_tail.ipol == 0.0) {
        it_end = 1;
    }

    for (Int_t i = 1; i <= it_end; ++i) {
        for (_tail.ita = 1; _tail.ita <= 2; ++_tail.ita) {
            std::cout << "********** ita: " << _tail.ita
                      << " *********" << std::endl;
            if (_tail.ita == 1) {
                Bornin();
                BorninTest(sibt);
                std::cout << "sib1" << fSib << std::endl;
                std::cout << "sibt" << sibt << std::endl;
                if (fSib == 0.0) {
                    fTai[1] = 0.;
                    continue;
                }
            }
            qqt(fTai[_tail.ita]);
            std::cout << "tai[" << _tail.ita
                      << "]\t"  << fTai[_tail.ita] << std::endl;
        }
         fDeltaInf = 0.;
         Double_t extai1 = TMath::Exp(kAlpha / kPi * fDeltaInf);
         fSig = fSib * extai1 * (1. + kAlpha / kPi * (fDelta - fDeltaInf)) +
                    fTai[1] + fTai[2];
    }

}



void TRadCor::Deltas(const Double_t massLepton2)
{
    Double_t sum = VacPol();

    Double_t xxh = _Sxy.s - _Sxy.y - _phi.vv10;
    Double_t ssh = _Sxy.x + _Sxy.y - _phi.vv20;

#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "sum    " << std::setw(20) << std::setprecision(10) << sum  << std::endl;
    std::cout << "xxh    " << std::setw(20) << std::setprecision(10) << xxh  << std::endl;
    std::cout << "ssh    " << std::setw(20) << std::setprecision(10) << ssh  << std::endl;
#endif

    Double_t alss = TMath::Power(ssh, 2) - 2. * _phi.p22 * (2 * massLepton2);
    Double_t alxx = TMath::Power(xxh, 2) - 2. * _phi.p22 * (2 * massLepton2);

    if (alss < 0) std::cout << "deltas: alss < 0 " << alss << std::endl;
    if (alxx < 0) std::cout << "deltas: alxx < 0 " << alxx << std::endl;

//    Double_t sqlss = TMath::Sqrt(TMath::Max(0., alss));
//    Double_t sqlxx = TMath::Sqrt(TMath::Max(0., alxx));

//    Double_t allss = TMath::Log((sqlss + ssh) / (-sqlss + ssh)) / sqlss ;
//    Double_t allxx = TMath::Log((sqlxx + xxh) / (-sqlxx + xxh)) / sqlxx ;

    Double_t dlm = TMath::Log(_Sxy.y / massLepton2);

//    Double_t sfpr = dlm * dlm / 2. - dlm * TMath::Log(ssh * xxh / massLepton2 / _phi.p22) -
//                    TMath::Power((TMath::Log(ssh / xxh)), 2) / 2. +
//                    HapradUtils::fspen(1 - _phi.p22 * _Sxy.y / ssh / xxh) - kPi * kPi / 3.;

    fDeltaInf = (dlm - 1) * TMath::Log(TMath::Power(_phi.p22 - kMassC2, 2) / ssh / xxh);
    fDelta    = fDeltaInf + sum +
                (1.5 * dlm - 2. - 0.5 * TMath::Power(TMath::Log(xxh / ssh),2) +
                        HapradUtils::fspen(1. - _phi.p22 * _Sxy.y / ssh / xxh) - kPi * kPi / 6.);
}



Double_t TRadCor::VacPol(void)
{
    Double_t leptonMass[3] = { 0.26110  * TMath::Power(10,-6),
                               0.111637 * TMath::Power(10,-1),
                               3.18301 };

    Double_t suml = 0;
    for (Int_t i = 0; i < 3; ++i) {
        Double_t a2    = 2 * leptonMass[i];
        Double_t sqlmi = TMath::Sqrt(_Sxy.y * _Sxy.y + 2 * a2 * _Sxy.y);
        Double_t allmi = TMath::Log((sqlmi + _Sxy.y) / (sqlmi - _Sxy.y)) / sqlmi;

        suml = suml + 2. * (_Sxy.y + a2) * allmi / 3. - 10. / 9. +
                    4. * a2 * (1. - a2 * allmi) / 3. / _Sxy.y;
    }

    Double_t a, b, c;

    if (_Sxy.y < 1) {
        a = -1.345 * TMath::Power(10,-9);
        b = -2.302 * TMath::Power(10,-3);
        c = 4.091;
    } else if (_Sxy.y < 64) {
        a = -1.512 * TMath::Power(10,-3);
        b = -2.822 * TMath::Power(10,-3);
        c =  1.218;
    } else {
        a = -1.1344 * TMath::Power(10,-3);
        b = -3.0680 * TMath::Power(10,-3);
        c =  9.9992 * TMath::Power(10,-1);
    }

    Double_t sumh;
    sumh = - (a + b * TMath::Log(1. + c * _Sxy.y)) * 2 * kPi / kAlpha;

#ifdef DEBUG
    std::cout << std::endl;
    std::cout.setf(std::ios::fixed);
    std::cout << "suml   " << std::setw(20) << std::setprecision(10) << suml  << std::endl;
    std::cout << "sumh   " << std::setw(20) << std::setprecision(10) << sumh  << std::endl;
    std::cout << std::endl;
#endif

    return suml + sumh;
}



void TRadCor::Bornin(void)
{

}



void TRadCor::BorninTest(Double_t& sib)
{

}



void TRadCor::qqt(Double_t& tai)
{

}
