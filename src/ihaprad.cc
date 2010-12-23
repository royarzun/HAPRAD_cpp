#include "haprad_constants.h"
#include "sxy.inc"
#include "tail.inc"
#include "phi.inc"
#include "epsil.inc"
#include "epsmarch.inc"
#include "TMath.h"
#include <algorithm>
#include <cstdlib>

using namespace TMath;

void ihaprad(Double_t& bmom, Int_t& ilepm, Int_t& iphi_radm, Double_t& epsphirm, Double_t& epstaum, Double_t& epsrrm, Double_t& xmas, Double_t& ymas, Double_t& zmas, Double_t& tmas, Double_t& phimas, Double_t& sib, Double_t& sig, Double_t& delinf, Double_t& delta, Double_t& tai[3])
{
    Double_t epspl;
    ipol = 0; //ipol - 1 -- long; 2 -- tran; 0 -- unpol <-- type of polarization
    iphi_had = 0; //iphi_had - integration over phi_{had} (1) or not (0)
    iphi_rad = iphi_radm; //iphi_rad - integration over phi_{rad} (1) or approximation (0)
    ilep = ilepm; //ilep - registered lepton: 1 - electron, 2 - muon
    epsphir = epsphirm; //relative accuracies of integration over phi_r (eq.3)
    epstau = epstaum; //relative accuracies of integration over tau (eq.3)
    epsrr = epsrrm; //relative accuracies of integration over R (eq.3)

    isf1 = 1.; //DIFINED TAIL.H
    isf2 = isf20; //DIFINED TAIL.H & HAPRAD_CONST.H isf20=4
    isf3 = 1.; //DIFINED TAIL.H
    un = 1.; //DIFINED TAIL.H
    pl = 0.; //DIFINED TAIL.H
    pn = 0.; //DIFINED TAIL.H
    qn = 0.; //DIFINED TAIL.H

    Double_t tmom = 0.; //tmom - momentum per nucleon
    Double_t snuc = 2. * amp * bmom; // S=2k_{1}p

    xs = xmas; //DIFINED SXY.H
    zdif = zmas; //DIFINED PHI.H
    tdif = tmas; //DIFINED PHI.H

    if (ymas > 0.) {
        ys = ymas; //DIFINED SXY.H
        Q2 = snuc * xs * ys; //DIFINED SXY.H Q2=xyS

    } else {
        Q2 = -ymas; // Q2
        ys = Q2 / (snuc * xs); //DIFINED SXY.H Bjorken Q2 ???

    }

    Double_t yma = 1. / (1. + Power(amp, 2.) * xs / snuc); // Maximum Bjorken Q2 ???
    amc2 = Power((amp + amhh), 2.); //DIFINED HAPRAD_CONST.H
    Double_t ymi = (amc2 - Power(amp, 2.)) / (snuc * (1. - xs)); //minimum Bjorken Q2 ???

    if (ys > yma || ys < ymi || xs > 1. || xs < 0.) {
        std::cout << " Warning! Wrong kinematics!!!! skip the point!" << std::endl;
        std::cout << " ys= " << ys << std::endl;
        std::cout << " xs= " << xs << std::endl;
        return 0;

    }

    conkin(snuc); //THE FUNCTION IS DEFINED HERE
    ehad = anu * zdif; //DIFINED PHI.H energy of hadron E_{h}=Nu*Z_{h}
    Double_t sqnuq = Sqrt(Power(anu, 2.) + Q2); //sqrt(Q^{2}+Nu^{2})

    if (ehad < amhh) { // the energy of detected hadron can't be less that its mass
        std::cout << " Warning! Wrong kinematics!!!! skip the point!" << std::endl;
        std::cout << " ehad =" << ehad << std::endl;
        return 0;
    }

    pph = Sqrt(Power(ehad, 2.) - Power(amhh, 2.)); //DIFINED PHI.H detected hadron momentum

    if (tdif > 0.) {
        pth = tdif; //DIFINED PHI.H

        if (pph < pth) { //hadron'S momentum can't be less than transverse momentum
            std::cout << " Warning! Wrong kinematics!!!! skip the point!" << std::endl;
            std::cout << " pph =" << pph << std::endl;
            std::cout << " pth =" << pth << std::endl;
            return 0;
        }
    }
    plh = Sqrt(Power(pph, 2.) - Power(pth, 2.)); //DIFINED PHI.H detected hadron'S longitudinal momentum

    if (pph > pth) an = an * sqrt_lq / 2. / amp / plh; //an DEFINED SXY.H
    else an = 0.;

    tdif = Power(amhh, 2.) - Q2 + 2. * (sqnuq * plh - anu * ehad);
    std::cout << "pl: " << plh << "\t" << tdif << "\t" << (plh - (tdif + Q2 - Power(amhh, 2.) + 2. * anu * ehad) / 2. / sqnuq) << std::endl;
    sibort = semi;
    return;
}

void deltas(Double_t &delta, Double_t &delinf, Double_t &tr)
{
//     sibor is born cross section with polarized initial
//     lepton and polarized target
//     siamm is contribution of anomalous magnetic moment.
    Double_t sqlxx;
    Double_t allss;
    Double_t allxx;
    Double_t dlm;
    Double_t sfpr;
    Double_t delta0;
    Double_t delta_old;
    Double_t aj0;
    Double_t xxh = S - Q2 - vv10;
    Double_t alss = Power(ssh, 2.) - 2. * p22 * al2;
    Double_t alxx = Power(xxh, 2.) - 2. * p22 * al2;

    if (alss < 0.) std::cout << "deltas: alss<0 " << alss << std::endl;

    Double_t sqlss = Sqrt(Max(0., alss));

    if (alxx < 0.) std::cout << "deltas: alxx<0" << alxx << std::endl;

    sqlxx = Sqrt(Max(0., alxx));
    allss = Log((sqlss + ssh) / (-sqlss + ssh)) / sqlss ;
    allxx = Log((sqlxx + xxh) / (-sqlxx + xxh)) / sqlxx ;
    dlm = Log(Q2 / aml2);
    sfpr = Power(dlm, 2.) / 2. - dlm * Log(ssh * xxh / aml2 / p22) - Power((Log(ssh / xxh)), 2.) / 2. + fspen(1. - p22 * Q2 / ssh / xxh) - pi2 / 3.;
    delta0 = (ssh * allss + xxh * allxx) / 2. + sfpr;
    delta_old = deltai + delta0 + del1 + del2 + sum;
    delinf = (dlm - 1.) * Power(Log(Power((p22 - amc2), 2.) / ssh / xxh));
    tr = alfa / pi * (dlm - 1.);
    delta = delinf + sum + (1.5 * dlm - 2. - 0.5 * Power(Log(xxh / ssh), 2.) + fspen(1. - p22 * Q2 / ssh / xxh) - pi2 / 6.);
    return;
}



void conkin(Double_t  snuc)
{
    //Set kinematics constants
    amt = amp;
    aml2 = amlep(ilep);
    aml2 = Power(amlep(ilep), 2.);
    al2 = 2. * aml2;
    pi2 = Power(Pi(), 2.);
    amh = amp;
    ap = 2. * amp;
    amp2 = amp * amp;
    ap2 = 2. * amp2;

    S = snuc * amp / amh;
    X = S * (1. - ys);
    S_x = S - X;
    S_p = S + X;
    ym = Q2 + al2;
    tpl = Power(S, 2.) + Power(X, 2.);
    tmi = Power(S, 2.) - Power(X, 2.);
    w2 = amp2 + S - Q2 - X;
    lambda_s = S * S - al2 * ap2;
    lambda_m = Q2 * Q2 + 4. * aml2 * Q2;
    lambda_q = Power(S_x, 2.) + 4. * amp2 * Q2;

    if (lambda_s < 0.) std::cout << "conkin: lambda_s < 0 " << lambda_s << std::endl;
    sqrt_ls = Sqrt(Max(0., lambda_s));
    if (lambda_x < 0.) std::cout << "conkin: lambda_x < 0 " << lambda_x << std::endl;
    sqrt_lx = Sqrt(Max(0., lambda_x));
    if (lambda_q < 0.) std::cout << "conkin: lambda_q < 0 " << lambda_q << std::endl;
    sqrt_lq = Sqrt(Max(0., lambda_q));
    if (lambda_m < 0.) std::cout << "conkin: lambda_m < 0 " << lambda_m << std::endl;
    sqlm = Sqrt(Max(0., lambda_m));

    allm = Log((sqlm + Q2) / (sqlm - Q2)) / sqlm;
    anu = S_x / ap;
    axy = Pi() * (S - X);
    an = Power(Pi(), 2.) * ys * S_x * amp / 2. / sqrt_lq * barn;

    taMax = (S_x + sqrt_lq) / ap2;
    tamin = -Q2 / amp2 / taMax;
    as = S / 2. / aml / sqrt_ls;
    bs = 0.;
    cs = -aml / sqrt_ls;

    if (ipol != 2.) {
        ae = amp / sqrt_ls;
        be = 0.;
        ce = -S / ap / sqrt_ls;
    } else {
        Double_t sqn_cmp = S * X * Q2 - lambda_q * aml2 - amp2 * Q2 * Q2;
        if (sqn_cmp > 0.) {
            sqn = Sqrt(sqn_cmp);
        } else {
            std::cout << "conkin: sqn=NaN " << sqn_cmp << std::endl;
            sqn = 0.;
        }
        ae = (-S * X + ap2 * ym) / sqrt_ls / sqn / 2.;
        be = sqrt_ls / sqn / 2.;
        ce = -(S * Q2 + al2 * S_x) / sqrt_ls / sqn / 2.;
    }
    apq = -Q2 * (ae - be) + ce * S_x;
    apn = (Q2 + 4. * aml2) * (ae + be) + ce * S_p;
    dk2ks = as * ym + al2 * bs + cs * X;
    dksp1 = as * S + bs * X + cs * ap2;
    dapks = 2. * (al2 * (as * ae + bs * be) + ap2 * cs * ce + ym * (as * be + bs * ae) + S * (as * ce + cs * ae) + X * (bs * ce + cs * be));
    return;
}



void bornin(Double_t& sibor)
{
    Double_t sfm0[8];
    Double_t tm[8];
    Double_t hi2;

    strf(0., 0., 0., sfm0); //llamada a strf
    hi2 = 0.25;

    tm[0] = Q2;
    tm[1] = (S * X - amp2 * Q2) / 2.;
    tm[2] = (vv10 * vv20 - Power(amhh, 2.) * Q2) / 2.;
    tm[3] = (S * vv20 + X * vv10 - zdif * S_x * Q2) / 2.;

    Double_t aa = S_x * (zdif - 2. * amp * plh / sqrt_lq) / 2. / amp2;
    Double_t ssum = 0.;
    for (Int_t i = isf1, i <= isf2, i = i + isf3) {
        ssum = ssum + tm[i] * sfm0[i];
        sibor = ssum * an / Q2 / Q2 * 2.;
    }
    return;
}


Double_t vacpol(Double_t t)
{
    //  am2 : squared masses of charge leptons
    Double_t am2[3] = {0.26110 * Power(10, -6), 0.111637 * Power(10, -1), 3.18301};
    Double_t a;
    Double_t sqlmi;
    Double_t aaa;
    Double_t bbb;
    Double_t ccc;
    Double_t sumh;
    Double_t suml = 0.;

    for (Int_t i = 0; i < 3; i++) {
        a2 = 2. * am2[i];
        sqlmi = Sqrt(t * t + 2. * a2 * t);
        allmi = Log((sqlmi + t) / (sqlmi - t)) / sqlmi;
        suml =  suml + 2. * (t + a2) * allmi / 3. - 10. / 9. + 4. * a2 * (1. - a2 * allmi) / 3. / t;
    }

    if (t < 1) {
        aaa = -1.345 * Power(10, -9);
        bbb = -2.302 * Power(10, -3);
        ccc = 4.091;
    } else if (t < 64) {
        aaa = -1.512 * Power(10, -3);
        bbb = -2.822 * Power(10, -3);
        ccc =  1.218;
    } else {
        aaa = -1.1344 * Power(10, -3);
        bbb = -3.0680 * Power(10, -3);
        ccc = 9.9992 * Power(10, -1);
    }
    sumh = -(aaa + bbb * Log(1. + ccc * t)) * 2. * pi / alfa;

    return suml + sumh;
}


void qqt(Double_t& tai)
{
    Double_t phiar[4];
    Double_t tar[6];
    Double_t ta1, ta2;
    Double_t ot, otr, am[2], bm[2], wrk[500], re, re2;
    extern Double_t qqtphi;

    if (ita == 1) {
        if (iphi_rad == 1) {
            simpsx(0, (2. * pi), 150, epsphir, qqtphi, tai);
            tai = an * alfa * tai / Power(pi, 2) / 4. / sqrt_lq;

        else if(iphi_rad == 0)
            tai = an * alfa / pi * qqtphi(0) / 2 / sqrt_lq;

        }

    } else {
        ta1 = -Q2 / S;
        ta2 = Q2 / X;
        phiar[0] = 0.;
        phiar[1] = 0.01 * pi;
        phiar[2] = 2. * pi - 0.01 * pi;
        phiar[3] = 2. * pi;
        tar[0] = tamin;
        tar[1] = ta1 - 0.15 * (ta1 - tamin);
        tar[2] = ta1 + 0.15 * (ta2 - ta1);
        tar[3] = ta2 - 0.15 * (ta2 - ta1);
        tar[4] = ta2 + 0.15 * (taMax - ta2);
        tar[5] = taMax;
        ot = Power(10, -3);

        Int_t id = 1;
        Double_t rere = 0.;
        for (Int_t iph = 0; iph < 4; iph++) {
            for (Int_t ico = 0; ico < 6; ico++) {
                am[0] = tar[ico];
                bm[0] = tar[ico+1];
                am[1] = phiar[iph];
                bm[1] = phiar[iph+1];

                if (am[1] > bm[1]) std::cout << am[1] << " < " << bm[1] << std::endl;
                if (am[0] > bm[0]) std::cout << am[0] << " < " << bm[1] << std::endl;

                Int_t mir = 10000;
                Int_t ma = 10 * mir;
                //      Integracion de d01fce
               // d01fce(2, am, bm, mir, ma, rv2tr, ot, otr, 500, wrk, re, id);
                write(*, '(1x,''tai:'',2i3,g15.6,f8.4,i9,i3)')ico, iph, re, otr, mir, id; // --->PENDIENTEEE
                rere = rere + re;
            }
        }
        tai = -alfa / (64. * Power(pi, 5.) * sqrt_lq * amp) * an * rere;
    }
    return;
}


Double_t rv2tr(Int_t ndim, Double_t *arg)
{

    Double_t tm[8][6], sfm[5], pres;
    Double_t vv;
    Double_t dmu;
    Double_t ta;
    if (ndim < 2 || ndim > 15) return;

    phirad = arg[1];
    ta = arg[0];
    amh2 = Power(amhh, 2);
    amu2 = Power(amhu, 2);
    vv = (1 - zdif) * S_x + tdif + amp2 - amu2;
    dmu = (ehad - plh * (S_x - ta * ap2) / sqrt_lq - ap2 * pth * Sqrt((ta - tamin) * (taMax - ta)) * Cos(phidif - phirad) / sqrt_lq) / amp;
    Double_t fwiw = 1. + ta - dmu;
    Double_t r = vv / fwiw;
    Double_t tldq2 = Q2 + r * ta;
    Double_t tldw2 = w2 - r * (1 + ta);
    Double_t tldtm = tdif - r * (ta - dmu);
    sffun(sfm, tldq2, tldw2, tldtm);
    tails(ta, tm, dmu);
    Double_t podinlz = 0;
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 4; j++) {
            pres = sfm[i] * tm[i][j] / Power(tldq2, 2) * Power(r, (j - 1)); // REVISAR POR EL J-1 ERA J-2 SEGUN LOGICA DE ARREGLOS ANTERIOR
            podinlz = podinlz + pres;
        }
    }
    rv2tr = podinlz / fwiw;
    return rv2tr;
}

void sffun(Double_t *sfm, Double_t q2, Double_t w2, Double_t t)
{
    //sfm[5]
    Double_t sffun;
    Double_t st, sl, stt, slt, sltp, sfm10, coetr;
    Double_t sqw2 = Sqrt(w2);
    Double_t S_x = w2 + q2 - amp2;
    Double_t lambda_q = Power(S_x, 2) + 4. * amp2 * q2;
    Double_t sxt = S_x + t + amp2 - amu2;
    Double_t tq = t + q2 - amh2;
    Double_t sffun_cmp = Power((W2 - amu2 - amh2), 2) - 4 * amu2 * amh2;

    if (sffun_cmp < 0) {
        std::cout << "sffun: sqlw=NaN " << sffun_cmp << std::endl;
    }

    Double_t sqlw = Sqrt(Max(0, sffun_cmp));
    Double_t ssffun_cmp = q2 * Power(sxt, 2) - sxt * S_x * tq - amp2 * Power(tq, 2) - amh2 * lambda_q;

    if (ssffun_cmp < 0) {
        std::cout << "ssffun: qll=NaN " << ssffun_cmp << std::endl;
    }

    Double_t sqll = Sqrt(Max(0,ssffun_cmp));
    Double_t cspion = (2 * tq * w2 + (S_x - 2 * q2) * (w2 + amh2 - amu2)) / sqlw / Sqrt(lambda_q);
//c Exclusive peak model (cross sections sigma_L,T,LT... from MAID2003)
    exclusive_model(q2, sqw2, cspion, st, sl, stt, slt, sltp);
    Double_t sfm20;
    Double_t sfm2tl;
    Double_t sfm4tl;
    Double_t sfm5tl;
    Double_t sfm4tt;
    Double_t sfm3tt;
    Double_t sfm2tt;
    Double_t sfm5tl;
// Structure functions

    if (lambda_q > 0 && sqll > 0 && sqlw > 0) {
        sfm10 = st - stt;
        sfm20 = 4. * (st + sl) * q2 / lambda_q;
        sfm2tl = 2. * slt * Sqrt(q2) * (-S_x * tq + 2. * q2 * sxt) / (lambda_q * sqll);
        sfm4tl = -slt * Sqrt(q2) / sqll;
        sfm4tt = -2. * stt * (-S_x * tq + 2. * q2 * sxt) / Power(sqll, 2.);
        sfm3tt = 2. * stt * lambda_q / Power(sqll, 2.);
        sfm2tt = 2. * stt * (Power((-S_x * tq + 2. * q2 * sxt), 2.) - 2. * q2 * Power(sqll, 2.)) / (lambda_q * Power(sqll, 2.));
        sfm5tl = -sltp * Sqrt(q2) / sqll;

        coetr = 16. * pi * (w2 - amp2) * w2 / (alfa * sqlw) / barn * 1000.;
        sfm[0] = coetr * sfm10;
        sfm[1] = coetr * (sfm20 + sfm2tl + sfm2tt);
        sfm[2] = coetr * sfm3tt;
        sfm[3] = coetr * (sfm4tl + sfm4tt);
    } else {
        sfm[0] = 0;
        sfm[1] = 0;
        sfm[2] = 0;
        sfm[3] = 0;
    }

    if (sfm[2] != sfm[2]) {
        std::cout << "sffun: " << coetr, st, sl, stt, slt, sltp, q2, lambda_q, S_x, sqw2, cspion, sxt, sqll << std::endl

         }
         return;
}

Double_t qqtphi(Double_t phi)
{
    //se borro la varable tlm[4] por no ser utiñ para nada
    Double_t tar[6], re;
    extern Double_t rv2ln;
    phirad = phi;
    Double_t ep = Power(1, -12);
    Double_t ta1 = -Q2 / S;
    Double_t ta2 = Q2 / X;
    tar[0] = tamin;
    tar[1] = ta1 - 0.15 * (ta1 - tamin);
    tar[2] = ta1 + 0.15 * (ta2 - ta1);
    tar[3] = ta2 - 0.15 * (ta2 - ta1);
    tar[4] = ta2 + 0.15 * (taMax - ta2);
    tar[5] = taMax;
    Double_t res = 0;
    for (Int_t i = 0; i < 6; i++) {
        simptx(Log(xs + tar[i]) + ep, Log(xs + tar[i+1]) - ep, 100, epstau, rv2ln, re));
        res = res + re;
    }
    return res;
}
/*****************************************************************************************************
void tails() Ocupa como variable global a "ita" Q2 verifica su estado

******************************************************************************************************/
void tails(Double_t ta, Double_t *tm[][6], Double_t amu)
{
    Double_t phka, phkb;
    Double_t bb;
    Double_t b1;
    Double_t b2;
    Double_t c1;
    Double_t c2;
    Double_t sqrtmb;
    Double_t z1;
    Double_t z2;
    Double_t bi12;
    Double_t bi1p12;
    Double_t bis;
    Double_t bir;
    Double_t b1i;
    Double_t b11i;

    if (iphi_rad == 0 && ita == 1) {
        b2 = (-lambda_q * ta + S_p * S_x * ta + 2. * S_p * Q2) / 2.;
        b1 = (-lambda_q * ta - S_p * S_x * ta - 2. * S_p * Q2) / 2.;
        c1 = -(4. * (amp2 * Power(ta, 2) - S_x * ta - Q2) * aml2 - Power((S * ta + Q2), 2));
        c2 = -(4. * (amp2 * Power(ta, 2) - S_x * ta - Q2) * aml2 - Power((ta * X - Q2), 2));
        bb = 1;

        if (c1 < 0) {
            std::cout << "tails: sc1=NaN " << c1 << std::endl;
        }

        Double_t sc1 = Sqrt(Max(0, c1));

        if (c2 < 0) {
            std::cout << "tails: sc2=NaN " << c2 << std::endl;

        }
        Double_t sc2 = Sqrt(Max(0, c2));
        bi12 = sqrt_lq * (S_p * (S_x * ta + 2. * Q2)) / (sc1 * sc2 * (sc1 + sc2));
        bi1pi2 = sqrt_lq / sc2 + sqrt_lq / sc1;
        bis = sqrt_lq * (-b1 / sc1 / c1 + b2 / sc2 / c2) * aml2;
        bir = sqrt_lq * (b2 / sc2 / c2 + b1 / sc1 / c1) * aml2;
        b1i = -sqrt_lq * b1 / lambda_q / sqrt_lq;
        b11i = sqrt_lq * (3. * Power(b1, 2) - lambda_q * c1) / 2. / Power(lambda_q, 2) / sqrt_lq;
    } else {
        Double_t sqrtmb_comp = (ta - tamin)*(taMax - ta)*(S * X * Q2 - Power(Q2, 2) * amp2 - aml2 * lambda_q);
        if ( sqrtmb_comp > 0) {
            sqrtmb = Sqrt(sqrtmb_comp);
        } else {
            std::cout << "tails: sqrtmb=NaN " << sqrtmb_comp << std::endl;
            sqrtmb = 0;
        }
        z1 = (Q2 * S_p + ta * (S * S_x + ap2 * Q2) - ap * Cos(phirad) * sqrtmb) / lambda_q;
        z2 = (Q2 * S_p + ta * (X * S_x - ap2 * Q2) - ap * Cos(phirad) * sqrtmb) / lambda_q;
        bb = 1. / sqrt_lq / pi;
        bi12 = bb / (z1 * z2);
        bi1pi2 = bb / z2 + bb / z1;
        bis = (bb / Power(z2, 2) + bb / Power(z1, 2)) * aml2;
        bir = (bb / Power(z2, 2) - bb / Power(z1, 2)) * aml2;
        b1i = bb * z1;
        b11i = bb * Power(z1, 2);
    }
    Double_t hi2 = bis - ym * bi12;
    Double_t amh2 = Power(amhh, 2);
    Double_t zh = zdif;
    Double_t vvp = (vv10 + vv20) / 2.;
    Double_t vvm = (vv10 - vv20) / 2.;
    tm[0][0] = 4. * Q2 * hi2;
    tm[0][1] = 4. * ta * hi2;
    tm[0][2] = -2. * (bi12 * Power(ta, 2) + 2. * bb);
    tm[1][0] = 2. * hi2 * (S * X - amp2 * Q2);
    tm[1][1] = (bi0pi2 * S_x * S_p - bi12 * ta * Power(S_p, 2) + 2. * bir * S_p + 2. * hi2 * (S_x - 2. * amp2 * ta)) / 2.;
    tm[1][2] = (bi02 * ta * (2. * amp2 * ta - S_x) - bi1pi2 * S_p + 4. * amp2 * bb) / 2.;
    tm[2][0] = 1. * hi2 * (vv10 * vv20 - amh2 * Q2);
    tm[2][1] = -2. * ((amh2 * ta - amu * vvm) * hi2 - bir * amu * vvp - bi0pi2 * vvm * vvp + bi12 * ta * Power(vvp, 2));
    tm[2][2] = bi01 * ta * (amh2 * ta - amu * vvm) - bi1pi2 * amu * vvp + 2. * amh2 * bb;
    tm[3][0] = -1. * (vv10 * S_x - 2. * S * vvp + Q2 * S_x * zh) * hi2;
    tm[3][1] = -2. * bi02 * ta * vvp * S_p + bi1pi2 * (S_p * vvm + S_x * vvp) + bir * (amu * S_p + 2. * vvp) + hi2 * (S_x * (amu - 2. * ta * zh) + 2. * vvm);
    tm[3][2] = bi01 * ta * (S_x * (ta * zh - amu / 2d0) - vvm) - bi1pi2 * (amu / 2. * S_p + vvp) + 2. * S_x * zh * bb;
    return;
}


Double_t rv2ln(Double_t taln, Bool_t nonzero)
{
    Double_t d2kvir;
    Double_t phka;
    Double_t phkb;
    Double_t factor;
    Double_t rmin;
    Double_t rMax;
    Double_t res, rv, DSIMPS, aval, bval;
//  Double_t vyv(0: 1000)// comentado dado que al parecer no se ocupa
    Bool_t nonzero;
    extern Double_t podinl;
    Double_t ta = Exp(taln) - Q2 / S_x;
    taa = ta;
    Double_t costk = (S_x - ap2 * ta) / sqrt_lq;
    Double_t sintk;

    if (abs(costk) <= 1) {
        sintk = Sqrt(1. - Power(costk, 2));

    } else {
        std::cout << "rv2ln: costk>1 " << costk << std::endl;
        sintk = 0;
    
	}

    d2kvir = (ehad - plh * costk - pth * sintk * Cos(phirad - phidif)) / amp;
    phka = 0.5 * (ehad - plh * costk) / amp;
    phkb = 0.5 * (-pth * sintk) / amp;
    daa = d2kvir;
    factor = 1. + ta - d2kvir;
    tails(ta, tm, d2kvir);
    if (ita == 1) {
        strf(0, 0, 0, sfm0);
        rmin = Power(10, -8);
        rMax = (p22 - amc2) / factor;
        if (nonzero) {
            // Check minimum r
            nonzero = false;
            for (Int_t i = 0; i < 5; i++) {
                if (sfm0[i] != 0) nonzero = true;
            }
            if (nonzero) {
                if ((ap * d2kvir - 2. * ehad) > 0) {
                    aval = (2. * ehad * S_x - ap * (Power(amhad(ivec), 2) - Q2 - tdif)) / (ap * d2kvir - 2. * ehad);
                    bval = Power(pph, 2) * 2. * Sqrt(Power(S_x, 2) + Power((ap * Q2), 2)) / (ap * d2kvir - 2 * ehad);
                    rmin = aval - bval;
                    rmin = (abs(aval) + abs(bval)) * Power(1, -7);

                } else {
                    aval = 2. * ehad * S_x - ap * (Power(amhad(ivec), 2) - Q2 - tdif);
                    bval = Power(pph, 2) * 2. * Sqrt(Power(S_x, 2) + Power((ap * Q2), 2));
                    rmin = (abs(aval) + abs(bval)) * Power(10, -7);

                }
            }
        }


    } else if (ita == 2) {
        res = podinl((p22 - amp2) / factor) / factor;

    }

    return res * (Q2 / S_x + ta);
}


Double_t  podinl(Double_t r)
{
//
//    integrand (over r )
    Double_t pp, pres, sfm[8], ta, d2kvir, rm
    Bool_t reget;
    ta = taa;
    d2kvir = daa;
    strf(ta, d2kvir, r, sfm);
    podinl = 0;

// Check lower bound
    if (ita == 1) {
        rm = r;
        while (1) {
            reget = false;
            for (Int_t i = 0; i < 5; i++) {
                if (sfm[i] == 0 && sfm0[i] != 0 && r < Power(10, -6)) {
                    reget = true;
                    std::cout <<  i, sfm[i], sfm0[i] << std::endl;
                }
            }
            if (reget) {
                rm = rm + Power(10, -8);
                std::cout << "regetting " << rm, r << std::endl;
                strf(ta, d2kvir, rm, sfm);
                /*stop*/
                continue;
            }
            break;
        }
    }

//      print*,'podinl start ',isf1,isf2,isf3,i1(isf)+i2(isf)-1

    for (Int_t isf = 0 ; isf <= isf2; isf++) {
        for (Int_t irr = 0; irr <= (i1[isf] + i2[isf] - 2); irr++) {
            pp = sfm[isf];
            if (irr == 1 && ita == 1) pp = pp - sfm0[isf] * Power((1 + r * ta / Q2), 2);

            pres = pp * Power(r, (irr - 2)) / Power((Q2 + r * ta), 2);
            podinl = podinl - tm[isf][irr] * pres;
//     if(irr.eq.1) print*,'pod: ',isf,irr,sfm0(isf),sfm(isf),pp,pres,r**(irr-2)
        }
    }
    return podinl;
}

void strf(Double_t ta, Double_t d2kvir, Double_t rr, Double_t  * sfm)
{
//c
//c     the programm calculates deep inelastic (ita=1),
//c     elastic (ita=2), quasielastic (ita=3) structure functions
//c     in kinematical point (ta,rr).
//c    rr=S_x-tt,
//c    ta=(t-Q2)/rr,
//c    where tt=t+amf2-amp2, amf2 is invarint mass of final hadrons
    Double_t tldq2;
    Double_t tldtd;
    Double_t tldaly;
    Double_t tldsqly;
    Double_t tldp22;
    Double_t phq;
    Double_t b1, b2, b3, b4, H1z, H2z, H3z, H4z, dum;
    Double_t tldnu;
    Double_t aks;
    Double_t zh;
    Double_t tldsq;
    Double_t tldplh;
    Double_t tldpt2;
    Double_t epsnu;
    Double_t epst;
    Double_t epsphq;
    Double_t epsq;
    Double_t epspl;
    Double_t epspt2;

    for (Int_t i = 0 ; i < 8 ; i++) {
        sfm[i] = 0;
    }

    tldq2 = Q2 + rr * ta;
    tldtd = tdif - rr * (ta - d2kvir);
    tldaly = Power((S_x - rr), 2) + 4. * amp2 * tldq2;
    tldsqly = Sqrt(tldaly);
    tldp22 = p22 - rr * (1. + ta - d2kvir);
    phq = (Power(amhh, 2) - tldq2 - tldtd) / 2.;

//      print*,'---------------Begin-----------------'
//      write(*,'(3f22.19)') S_x,rr,ap
    tldnu = (S_x - rr) / ap;
    aks = tldq2 / ap / tldnu;
    zh = ehad / tldnu;
    tldsq = Sqrt(tldq2 + Power(tldnu, 2));
//c      tldplh =(tldtd + tldq2 - Power(amhad(ivec),2) + 2 * tldnu * ehad) / 2 / tldsq;
    tldplh = (ehad * tldnu - phq) / tldsq;
    tldpt2 = Power(pph, 2) - Power(tldplh, 2);
    epsnu = Sqrt(Power((S_x * epsmarch / ap), 2) + Power((rr * epsmarch / ap), 2) + Power(((S_x - rr) / ap * epsmarch), 2));
    epst = Sqrt(Power((tdif * epsmarch), 2) + (Power(((ta - d2kvir) * rr * epsmarch), 2) + Power((rr * ta * epsmarch), 2) + Power((rr * d2kvir * epsmarch), 2)));
    epsphq = Sqrt(Power((2 * Power(amhh, 2) * epsmarch), 2) + Power((tldq2 * epsmarch), 2) + Power(epst, 2)) / 2;
    epsq = 1. / 2. / Sqrt(tldq2 + Power(tldnu, 2)) * Sqrt(Power((tldq2 * epsmarch), 2) + Power((2 * (tldnu, 2) * epsnu), 2));
    epspl = Sqrt(Power((ehad * epsmarch * tldnu / tldsq), 2) + Power((ehad * epsnu / tldsq), 2) + Power((epsphq / tldsq), 2) + Power(((ehad * tldnu - phq) / Power(tldsq, 2) * epsq), 2));
    epspt2 = 2. * Sqrt(Power(pph, 4) * Power(epsmarch, 2) + Power(tldplh, 4) * Power(epspl, 2));
//c      print*,'pt2= ',tldpt2,epspt2,(tldpt2**2-epspt2**2),epspl,tldplh,pph
//c      if(tldpt2.lt.0.d0) return
    if (tldpt2 < 0 && (Power(tldpt2, 2) - Power(epspt2, 2)) > 0) return;

//c       b1=0.d0
//c       b2=0.d0
//c       b3=0.d0
//c       b4=0.d0
//c print *,'q2',tldq2
//c Call semi-inclusive model (H_i defined in Mulders)

//c      dum=tldq2/ap/aks
//c      if(tldnu.ne.dum) then
//c      print*,'redefine'
//c      aks=tldq2/ap/tldnu
//c      zh=ehad/tldnu
//c      tldsq=sqrt(tldq2+tldnu**2)
//c      tldplh=(ehad*tldnu-phq)/tldsq
//c      tldpt2=pph**2-tldplh**2
//c      if(tldpt2.lt.0.d0) return
//c      write(*,'(10f22.19)') tldq2,aks,zh,tldpt,tldnu,Ehad,pph,
//c     &dum,(tldq2/ap/tldnu),(zh*tldnu)
//c      endif

//c      write(*,'(10f22.19)') tldq2,aks,zh,tldpt,tldnu,Ehad,pph,
//c     &dum,(tldq2/ap/tldnu),(zh*tldnu)

//c Recover Q2
    Double_t a = S / ap * (S / ap - anu) * tldq2 / Q2;
    Double_t tlde = (tldnu + Sqrt(Power(tldnu, 2) + 4 * a)) / 2;
    Double_t tldy = tldnu / tlde;

    semi_inclusive_model(tldq2, aks, tldy, zh, tldpt2, tldp22, tldplh, H1z, H2z, H3z, H4z);
//c       print*,'----------------------------------------------'

//c      print*,'Hs ',tldq2,aks,zh,tldpt,H1z,H2z,H3z,H4z

//c print*,'semi-inclusive'
//c print*,Eb,tldq2,aks,zh,sqrt(tldpt2)
//c print*,H1z,H2z,H3z,H4z


//c No photon emission (k-->0)
//c      aa=S_x*(zdif-2.*amp*plh/sqrt_lq)/2.d0/amp2
//c h(1)=  zdif*(2.d0*Q2*xs*sfm0(1)-pth**2*sfm0(4))/amp/Q2/xs
//c h(2)=2.*zdif*(2.*xs*(xs*sfm0(2)+aa*sfm0(3))*lambda_q
//c     . +sfm0(4)*(aa**2*lambda_q-2.*pth**2*Q2))/amp/xs/Q2/lambda_q
//c        h(3)=2.*sfm0(4)*zdif/amp/Q2/xs
//c h4=-2.*zdif*(xs*sfm0(3)+aa*sfm0(4))/amp/Q2/xs

//c Including kinematic shift due to photon emission
    Double_t epsi = ap2 / (S_x - rr);
    Double_t aa = (S_x - rr) * (zh - 2 * amp * tldplh / tldsqly) / 2. / amp2;
    Double_t h1 = zh * (2. * tldq2 * aks * h1z - tldpt2 * h4z) / amp / tldq2 / aks;
    Double_t h2 = 2 * zh * (2 * aks * (aks * h2z + aa * h3z) * tldaly + h4z * (Power(aa, 2) * tldaly - 2 * tldpt2 * tldq2)) / amp / aks / tldq2 / tldaly;
    Double_t h3 = 2 * h4z * zh / amp / tldq2 / aks;
    Double_t h4 = -2 * zh * (aks * h3z + aa * h4z) / amp / tldq2 / aks;

//c      print*,epsi,aa,h1,h2,h3,h4


    sfm[0] = un * h1;
    sfm[1] = un * h2;
    sfm[2] = un * h3;
    sfm[3] = un * h4;
//c      sfm(5)=epsi**2*b1
//c      sfm(6)=epsi**3*(b2/3.d0+b3+b4)
//c      sfm(7)=epsi*(b2/3.d0-b3)
//c      sfm(8)=epsi**2*(b2/3.d0-b4)


//c      print*,'strf: ',tldq2,aks,zh,tldpt2,tldpt,H1z,H2z,H3z,H4z,
//c     &tldnu,S_x,ap,ehAD

    return;
}
