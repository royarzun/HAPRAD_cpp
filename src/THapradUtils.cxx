#include "THapradUtils.h"
#include "TMath.h"

extern "C" {
    void init_pdf_(int&, int&);

    void exec_structm_(double&, double&, double&, double&, double&, double&,
                       double&, double&, double&, double&, double&);

    void exec_pkhff(int&, int&, double&, double&, double*, double*, double*,
                    double*, double*, double*);
}

namespace HapradUtils {

    Double_t fspen(const Double_t x)
    {
        const Double_t f1 = 1.644934;
        Double_t result;

        if (x < -1) {
            Double_t logprod;
            logprod = TMath::Log(1 - x) * TMath::Log(x * x / (1 - x));
            result = - 0.5 * logprod - f1 + fspens(1 / (1 - x));
        } else if (x < 0) {
            Double_t log2;
            log2 = TMath::Log(1 - x) * TMath::Log(1 - x);
            result = - 0.5 * log2 - fspens(x / (x - 1));
        } else if (x <= 0.5) {
            result = fspens(x);
        } else if (x <= 1) {
            Double_t logprod;
            logprod = TMath::Log(x) * TMath::Log(1 - x + TMath::Power(10, -10));
            result = f1 - logprod - fspens(1 - x);
        } else if (x <= 2) {
            Double_t logprod;
            logprod = TMath::Log(x) * TMath::Log((x - 1) * (x - 1) / x);
            result = f1 - 0.5  * logprod + fspens(1 - 1 / x);
        } else {
            Double_t log2;
            log2 = TMath::Log(x) * TMath::Log(x);
            result = 2 * f1 - 0.5 * log2 - fspens(1 / x);
        }

        return result;
    }



    Double_t fspens(const Double_t x)
    {
        Double_t f   = 0;
        Double_t a   = 1;
        Double_t an  = 0;
        Double_t tch = TMath::Power(10, -16);
        Double_t b;

        do {
            an = an + 1;
            a = a * x;
            b = a / TMath::Power(an, 2);
            f = f + b;
        } while (b - tch > 0);

        return f;
    }



    void SemiInclusiveModel(Double_t q2, Double_t X,
                            Double_t Y, Double_t Z,
                            Double_t pt2, Double_t mx2,
                            Double_t pl, Double_t& H1z,
                            Double_t& H2z, Double_t& H3z,
                            Double_t& H4z)
    {
        Double_t ac = 1.2025 * Power(10, -10);
        Double_t bc = -5.2703 * Power(10, -2);
        Double_t cc = 3.7467 * Power(10, -1);
        Double_t dc = 6.5397 * Power(10, -2);
        Double_t ec = -2.2136 * Power(10, -1);
        Double_t fc = -1.0621 * Power(10, -1);
        Double_t GTMD, SCALE, UPV, DNV, USEA, DSEA, STR, CHM, BOT, TOP, GD, XD, GL;
        Double_t r, xi;
        Int_t GPDF = 5;
        Int_t SPDF = 5;
        Int_t ISET = 1;
        Int_t ICHARGE = 1;
        Int_t nc = 0;
        /*Check the kinematics*/
        if (X < 0 || X > 1) return;

        if (Z < 0 || Z > 1) return;

        nc++;
        r = Sqrt(1. + Power(2 * mp * X, 2) / q2);
        xi = 2.*X / (1 + r);

        if (q2 > 1) SCALE = Sqrt(q2);
        else SCALE = 1.0000001;

        if (nc == 1) init_pdf_(GPDF, SPDF);

        GL = 1.033;

        exec_structm_(XD, SCALE, UPV, DNV, USEA, DSEA, STR, CHM, BOT, TOP, GL);

        Double_t ZD;
        if (Z > 0.01) ZD = Z;
        else ZD = 0.01000001;

        Double_t Q2D;
        if (q2 > 1) Q2D = q2;
        else Q2D = 1.0000001;

        Int_t IFINI;
        if (nc == 1) IFINI = 0;
        else IFINI = 1;

        Double_t uff[2], dff[2], sff[2], cff[2], bff[2], gff[2];

        exec_pkhff_(ISET, ICHARGE, ZD, Q2D, uff, dff, sff, cff, bff, gff);

        Double_t sgmpt = ac + bc * X + cc * Z + dc * Power(Z, 2) + ec * Power(Z, 2) + fc * X * Z;

        if (sgmpt < 0.02) sgmpt = 0.02;
        if (sgmpt > 0.15) sgmpt = 0.15;

        Double_t GMTD;

        if (pl > 0.15)
            GMTD  = Exp(-pt2 / (2.*sgmpt)) / (2.*pi * sgmpt);
        else
            GTMD  = Exp(-(pt2 + 2.*Power(pl, 2)) / (2.*sgmpt)) / (2.* pi * sgmpt);

        if (mx2 < Power((mp + mpi), 2)) return;

        Double_t uq, dq, sq, cq, bq, tq, gg;
        Double_t pi_thresh;
        Double_t H1, H2, H3m, H4m;
        Double_t xv, zv, q2v;
        Double_t rlt = 0.14;

        pi_thresh = Sqrt(1. - Power((mp + mpi), 2) / mx2);
        uq = eu2 * ((UPV + USEA) * uff[0] + USEA * uff[1]);
        dq = ed2 * ((DNV + DSEA) * dff[0] + DSEA * dff[1]);
        sq = es2 * (STR * sff[0] + STR * sff[1]);
        cq = ec2 * (CHM * cff[0] + CHM * cff[1]);
        bq = eb2 * (BOT * bff[0] + BOT * bff[1]);
        tq = 0.0;
        gg = 0.0;
        H2 = (uq + dq + sq + cq + bq + tq + gg) * GTMD * pi_thresh;
        H1 = H2 / (2. * X * (1. + rlt)) * (1. + 4. * Power(mp, 2) * Power(X, 2) / q2);
        xv = X;
        zv = Z;
        q2v = q2;

        if (X < 0.1) xv = 0.1;
        if (q2 < 1.) q2v = 1.;
        if (Z < 0.1) zv = 0.1;

        Double_t Ebeam, rt, rtz, cterm, m_cos_phi, m_cos_2phi;

        H3m = h3(X, q2, Z) * GTMD * pi_thresh;
        H4m = h4(X, q2, Z) * GTMD * pi_thresh;

// Check that <cos(phi)> and <cos(2phi)> are < 1

        Ebeam = q2 / (2. * mp * X * Y);
        rt = 1. - Y - mp * X * Y / (2. * Ebeam);
        rtz = Sqrt(rt / (1. + 2. * mp * X / (Y * Ebeam)));
        cterm = X * Power(Y, 2) * H1 + rt * H2;
        m_cos_phi = Sqrt(pt2 / q2) * (2. - Y) * rtz * H3m / (2.*cterm);

        if (Abs(4. * m_cos_phi) > 0.9)
            H3m = 0.9 * Sign(1., m_cos_phi) * (2.*cterm) / (Sqrt(pt2 / q2) * (2. - Y) * rtz) / 4.;
        m_cos_2phi = pt2 / q2 * Power(rtz, 2) * H4m / (2.*cterm);

        if (Abs(4.*m_cos_2phi) > 0.9)
            H4m = 0.9 * Sign(1., m_cos_2phi) * (2.*cterm) / (pt2 / q2 * Power(rtz, 2)) / 4.;

        H1z = H1;
        H2z = H2;
        H3z = H3m;
        H4z = H4m;

        return;
    }


    void dfint(Int_t narg, Double_t *arg, Double_t *nent, Double_t *ent, Double_t *table)
    {
        Double_t kd = 1;
        Double_t m = 1;
        Int_t ncomb[narg];
        Int_t ja = 0;
        Int_t jb, k;
        Int_t i = 0;

        while (i <= narg) {
            ncomb[i] = 1;
            jb = ja + nent[i]; // VERIFICAR SI LOS LIMITES ESTAN CORRECTOS!!!!
            Int_t j = ja;
            while (j <= jb) {
                if (arg[i] <= ent[j]) {
                    if (j != ja) {
                        jr = j - 1;
                        d[i] = (ent[j] - arg[i]) / (ent[j] - ent[jr]);
                        ient[i] = j - ja;
                        kd = kd + ient[i] * m;
                        m = m * nent[i];
                        ja = jb + 1;
                        break;
                        
                    } else {
                        j++;
                        jr = j - 1;
                        d[i] = (ent[j] - arg[i]) / (ent[j] - ent[jr]);
                        ient[i] = j - ja;
                        kd = kd + ent[i] * m;
                        m = m * nent[i];
                        ja = jb + 1;
                        break;

                    }

                }
                j++;
            }
            j = jb;
            j++;
            d[i] = (ent[j] - arg[i]) / (ent[j] - ent[jr]);
            ient[i] = j - ja;
            kd = kd + ient[i] * m;
            m = m * nent[i];
            ja = jb + 1;
            i++;
        }
        dfint = 0.;
        while (1) {
            fac = 1.;
            iadr = kd;
            ifadr = 1.;
            //cambiado a 0 dado el indice de los arreglos en fortran
            i = 0;
            while (i <= narg) {
                if (ncomb[i] == 0) {
                    fac = fac * d[i];
                    iadr  = iadr - ifadr;
                    ifadr = ifadr * nent[i];

                } else {
                    fac = fac * (1. -  d[i]);
                    ifadr = ifadr * nent[i];
                }
                i++;
            }
            dfint = dfint + fac * table[adr];
            il = narg;

            while (ncomb[il] == 0) {
                il--;
                if (il == 0) return;
            }

            ncomb[il] = 0;
            if (il == narg) {
                continue;
            } else {
                il++;
                k = il;
                while (k <= narg) {
                    ncomb[k] = 1;
                    k++;
                }
                continue;
            }
        }

    }

}
