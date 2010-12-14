#include "TRadCor.h"
#include "haprad_constants.h"
#include <iostream>


TRadCor::TRadCor()
: fEbeam(0),
  fX(0), fY(0), fZ(0), fPt(0), fPhi(0),
  fMaxMx2(0), fMx2(0), fSib(0), fSig(0), fDelta(0), fTail(0),
  _eps(), _phi(), _pol(), _Sxy(), _tail()
{
    // Default constructor
}



TRadCor::TRadCor(Double_t Ebeam, Double_t x, Double_t q2, Double_t z,
                 Double_t pt, Double_t phi, Double_t maxMx2)
: fEbeam(Ebeam),
  fX(x), fY(-q2), fZ(z), fPt(pt), fPhi(phi/kRadianDeg),
  fMaxMx2(0), fMx2(0), fSib(0), fSig(0), fDelta(0), fTail(0),
  _eps(), _phi(), _pol(), _Sxy(), _tail()
{
    // Normal constructor for a radiative correction object
    //
    // The Ebeam parameter is the energy of the beam, the x, q2, z, pt and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    Setup();
}



TRadCor::~TRadCor()
{
    // Default destructor
}



void TRadCor::SetParameters(Double_t Ebeam, Double_t x, Double_t q2, Double_t z,
                            Double_t pt, Double_t phi, Double_t maxMx2)
{
    // Set the values for the variables used to calculate the radiative
    // correction.
    //
    // The Ebeam parameter is the energy of the beam, the x, q2, z, pt and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    fEbeam = Ebeam;
    fX = x;
    fY = -q2;
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



void TRadCor::SetQ2(Double_t q2)
{
    fY = -q2;
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
    // Calculate the missing mass

    Double_t nu;
    Double_t Sx;

    nu = - fY / ( 2.0 * kMassProton * fX);
    Sx = 2.0 * kMassProton * nu;

    fMx2 = TMath::Power(kMassProton,2) + Sx * (1.0 - fZ) + fPt;
}



Double_t TRadCor::GetRCFactor(void)
{
    // Get the radiative correction factor. You must set the parameters before
    // using this method.

    Haprad();
    fRCFac = fSig / fSib;
    return fRCFac;
}



Double_t TRadCor::GetRCFactor(Double_t Ebeam, Double_t x, Double_t q2, Double_t z,
                              Double_t pt, Double_t phi, Double_t maxMx2)
{
    // Get the radiative correction factor for the given parameters.
    //
    // The Ebeam parameter is the energy of the beam, the x, q2, z, pt and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    SetParameters(Ebeam,x,q2,z,pt,phi,maxMx2);
    return GetRCFactor();
}



void TRadCor::Haprad(void)
{
    _tail.isf1 = 1;
    _tail.isf2 = 4;
    _tail.isf3 = 1;

    _tail.un = 1.;
    _tail.pl = 0.;
    _tail.pn = 0.;
    _tail.qn = 0.;

    _Sxy.xs = fX;
    _phi.zdif = fZ;
    _phi.tdif = fPt;

    Double_t snuc = 2. * kMassProton * fEbeam;        // S = 2k_{1}p

    if (fY >= 0.) {
        _Sxy.ys = fY;
        _Sxy.y = snuc * _Sxy.xs * _Sxy.ys;
    } else {
        _Sxy.y = - fY;
        _Sxy.ys = _Sxy.y / (snuc * _Sxy.xs);
    }

    Double_t mp2 = TMath::Power(kMassProton, 2);
    Double_t yma = 1. / (1. + mp2 * _Sxy.xs / snuc);
    Double_t ymi = (kMassC2 - mp2) / (snuc * (1. - _Sxy.xs));

    if (_Sxy.ys > yma || _Sxy.ys < ymi || _Sxy.xs > 1. || _Sxy.xs < 0.) {
        std::cout << " Warning! Wrong kinematics!!!! skip the point!"
                  << std::endl
                  << " ys= " << _Sxy.ys << std::endl
                  << " xs= " << _Sxy.xs << std::endl;
        return;
    }

    Conkin(snuc);

    _phi.ehad = _Sxy.anu * _phi.zdif;
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
        _phi.pth = TMath::Sqrt(_phi.pph * _phi.pph - _phi.pph * _phi.pph);
    }

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



void TRadCor::Conkin(const Double_t snuc)
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

    _Sxy.s   = snuc;
    _Sxy.x   = _Sxy.s * (1 - _Sxy.ys);
    _Sxy.sx  = _Sxy.s - _Sxy.x;
    _Sxy.sxp = _Sxy.s + _Sxy.x;
    _Sxy.ym  = _Sxy.y + 2 * ml2;
    _Sxy.tpl = TMath::Power(_Sxy.s, 2) + TMath::Power(_Sxy.x, 2);
    _Sxy.tmi = TMath::Power(_Sxy.s, 2) - TMath::Power(_Sxy.x, 2);
    _Sxy.w2  = mp2 + _Sxy.s - _Sxy.y - _Sxy.x;
    _Sxy.als = _Sxy.s * _Sxy.s - 2 * ml2 * (2 * mp2);
    _Sxy.alm = _Sxy.y * _Sxy.y + 4 * ml2 * _Sxy.y;
    _Sxy.aly = TMath::Power(_Sxy.sx, 2) + 4 * mp2 * _Sxy.y;

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

    _Sxy.tamax = (_Sxy.sx + _Sxy.sqly) / (2 * mp2);
    _Sxy.tamin = - _Sxy.y / mp2 / _Sxy.tamax;

    _pol.as = _Sxy.s / 2 / massLepton / _Sxy.sqls;
    _pol.bs = 0;
    _pol.cs = - massLepton / _Sxy.sqls;

    if (_tail.ipol != 2) {
        _pol.ae = kMassProton / _Sxy.sqls;
        _pol.be = 0;
        _pol.ce = - _Sxy.s / (2 * kMassProton) / _Sxy.sqls;
    } else {
        Double_t sqn;

        sqn = _Sxy.s * _Sxy.x * _Sxy.y - _Sxy.aly * ml2 - mp2 * _Sxy.y * _Sxy.y;
        if (sqn > 0) {
            sqn = TMath::Sqrt(sqn);
        } else {
            std::cout << "conkin: sqn = NaN " << sqn << std::endl;
            sqn = 0;
        }

        _pol.ae = (-_Sxy.s * _Sxy.x + (2 * mp2) * _Sxy.ym) / _Sxy.sqls / sqn / 2;
        _pol.be = _Sxy.sqls / sqn / 2;
        _pol.ce = -(_Sxy.s * _Sxy.y + (2 * ml2) * _Sxy.sx) / _Sxy.sqls / sqn / 2;
    }

    _pol.apq = -_Sxy.y * (_pol.ae - _pol.be) + _pol.ce * _Sxy.sx;
    _pol.apn = (_Sxy.y + 4 * ml2) * (_pol.ae + _pol.be) + _pol.ce * _Sxy.sxp;

    _pol.dk2ks = _pol.as * _Sxy.ym + (2 * ml2) * _pol.bs + _pol.cs * _Sxy.x;
    _pol.dksp1 = _pol.as * _Sxy.s + _pol.bs * _Sxy.x + _pol.cs * (2 * mp2);
    _pol.dapks = (2 * ml2) * (_pol.as * _pol.ae + _pol.bs * _pol.be) +
                 (2 * mp2) * _pol.cs * _pol.ce +
                 _Sxy.ym * (_pol.as * _pol.be + _pol.bs * _pol.ae) +
                 _Sxy.s  * (_pol.as * _pol.ce + _pol.cs * _pol.ae) +
                 _Sxy.x  * (_pol.bs * _pol.ce + _pol.cs * _pol.be);
    _pol.dapks *= 2.;

    return;
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

    Double_t costs, costx, sint;

    costx = (_Sxy.x * (_Sxy.s - _Sxy.x) - 2. * mp2 * _Sxy.y) / _Sxy.sqlx / _Sxy.sqly;
    costs = (_Sxy.s * (_Sxy.s - _Sxy.x) + 2. * mp2 * _Sxy.y) / _Sxy.sqls / _Sxy.sqly;

    Double_t sxy = _Sxy.s * _Sxy.x * _Sxy.y;
    Double_t y2  = _Sxy.y * _Sxy.y;

    sint = sxy - mp2 * y2 - ml2 * _Sxy.aly;

    if (sint > 0) {
        sint = 2. * kMassProton * TMath::Sqrt(sint) / _Sxy.sqls / _Sxy.sqly;
    } else {
        std::cout << "sphi: sint = NaN " << sint << std::endl;
        sint = 0.;
    }

    _phi.phk12 = - 0.5 * _Sxy.sqls * sint * _phi.pth / kMassProton;

    Double_t vv10 = costs * _phi.plh + sint * _phi.pth * TMath::Cos(fPhi);
    _phi.vv10 = (_Sxy.s * _phi.ehad - _Sxy.sqls * vv10) / kMassProton;

    Double_t vv20 = costx * _phi.plh + sint * _phi.pth * TMath::Cos(fPhi);
    _phi.vv20 = (_Sxy.s * _phi.ehad - _Sxy.sqls * vv20) / kMassProton;

    Double_t phk1, phk2;
    phk1 = 0.5 * (_Sxy.s * _phi.ehad - _Sxy.sqls * costs * _phi.plh) / kMassProton;
    phk2 = 0.5 * (_Sxy.s * _phi.ehad - _Sxy.sqlx * costx * _phi.plh) / kMassProton;

    _phi.phkp = phk1 + phk2;
    _phi.phkm = phk2 - phk1;

    Deltas();

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



void TRadCor::Deltas(void)
{

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
