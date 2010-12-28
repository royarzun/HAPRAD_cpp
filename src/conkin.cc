//CDECK  ID>, CONKIN.
//set up of kinematical observables
//****************** conkin *************************************
#include "haprad_consts.h"
#include "sxy.h"
#include "tail.h"
#include "phi.h"
#include "pol.h"

using namespace TMath;

void conkin(Double_t& snuc)
{
    // set of kinematical constants
    amt = amp; //DEFINED HAPRAD_CONST.H amp=0.938272
    aml = amlep[ilep-1]; //DEFINED HAPRAD_CONST.H amlep[2] = {0.511000e-3, 0.10565}
    aml2 = Power(amlep[ilep-1], 2); //DEFINED HAPRAD_CONST.H amlep[2] = {0.511000e-3, 0.10565}
    al2 = 2. * aml2; //DEFINED HAPRAD_CONST.H
    pi2 = Power(pi, 2); //DEFINED HAPRAD_CONST.H
    amh = amp; //DEFINED HAPRAD_CONST.H amp=0.938272
    ap = 2. * amp; //DEFINED HAPRAD_CONST.H
    amp2 = Power(amp, 2); //DEFINED HAPRAD_CONST.H
    ap2 = 2. * amp2; //DEFINED HAPRAD_CONST.H

    s = snuc * amp / amh; //DEFINED SXY.H amp / amh = 1
    x = s * (1. - ys); //DEFINED SXY.H X=(1-y)S
    sx = s - x; //DEFINED SXY.H S_{x}=S-X
    sxp = s + x; //DEFINED SXY.H S_{p}=S+X
    ym = y + al2; //DEFINED SXY.H
    tpl = Power(s, 2) + Power(x, 2); //DEFINED SXY.H
    tmi = Power(s, 2) - Power(x, 2); //DEFINED SXY.H
    w2 = amp2 + s - y - x; //DEFINED SXY.H
    als = s * s - al2 * ap2; //DEFINED SXY.H lambda_{s}=S^{2}-4m^{2}_{e}M^{2}_{p}
    alx = x * x - al2 * ap2; //DEFINED SXY.H lambda_{x}=X^{2}-4m^{2}_{e}M^{2}_{p}
    alm = y * y + 4. * aml2 * y; //DEFINED SXY.H
    aly = Power(sx, 2) + 4. * amp2 * y; //DEFINED SXY.H lambda_{q}=S^{2}_{x}+4M^{2}Q^{2}

    if (als < 0.) std::cout << "conkin: als<0 " << als << std::endl;
    sqls = Sqrt(Max(0., als)); //DEFINED SXY.H Sqrt(lambda_{s})
    if (alx < 0.) std::cout << "conkin: alx<0 " << alx << std::endl;
    sqlx = Sqrt(Max(0., alx)); //DEFINED SXY.H Sqrt(lambda_{x})
    if (aly < 0.) std::cout << "conkin: aly<0 " << aly << std::endl;
    sqly = Sqrt(Max(0., aly)); //DEFINED SXY.H Sqrt(lambda_{q})
    if (alm < 0.) std::cout << "conkin: alm<0 " << alm << std::endl;
    sqlm = Sqrt(Max(0., alm));  //DEFINED SXY.H
    allm = log((sqlm + y) / (sqlm - y)) / sqlm; //DEFINED SXY.H

    anu = sx / ap; //DEFINED SXY.H Nu=S_{x}/(2M_{p})

    Double_t axy = pi * (s - x); //pi*S_{x}
    //      an =2. * alfa**2 / sqls * axy * barn * amh / amp;
    an = pi * Power(alfa, 2) * ys * sx * amp / 2. / sqly * barn; //DEFINED SXY.H

    //      tamin = (sx - sqly) / ap2; //DEFINED SXY.H tau_{min}=(S_{x}+Sqrt(lambda_{q}))/2M^{2}_{p}
    taMax = (sx + sqly) / ap2; //DEFINED SXY.H tau_{max}=(S_{x}+Sqrt(lambda_{q}))/2M^{2}_{p}
    tamin = -y / amp2 / taMax; //DEFINED SXY.H tau_{min}
    as = s / 2. / aml / sqls; //DEFINED POL.H ERROR aml doesn't initiated ??????
    bs = 0.; //DEFINED POL.H
    cs = -aml / sqls; //DEFINED POL.H ERROR aml doesn't initiated ??????

    Double_t sqn;

    if (ipol != 2) {
        ae = amp / sqls; //DEFINED POL.H
        be = 0.; //DEFINED POL.H
        ce = -s / ap / sqls; //DEFINED POL.H

    } else {
        if ((s * x * y - aly * aml2 - amp2 * y * y) > 0.) {
            sqn = Sqrt(s * x * y - aly * aml2 - amp2 * y * y);

        } else {
            std::cout << "conkin: sqn=NaN " << (s * x * y - aly * aml2 - amp2 * y * y) << std::endl;
            sqn = 0.;

        }

        ae = (-s * x + ap2 * ym) / sqls / sqn / 2.;  //DEFINED POL.H
        be = sqls / sqn / 2.; //DEFINED POL.H
        ce = -(s * y + al2 * sx) / sqls / sqn / 2.;  //DEFINED POL.H
    }

    apq = -y * (ae - be) + ce * sx; //DEFINED POL.H
    apn = (y + 4. * aml2) * (ae + be) + ce * sxp; //DEFINED POL.H
    dk2ks = as * ym + al2 * bs + cs * x; //DEFINED POL.H
    dksp1 = as * s + bs * x + cs * ap2; //DEFINED POL.H
    dapks = 2. * (al2 * (as * ae + bs * be) + ap2 * cs * ce + ym * (as * be + bs * ae) + s * (as * ce + cs * ae) + x * (bs * ce + cs * be)); //DEFINED POL.H
    return 0;
}
