#include "THapradUtils.h"
#include "haprad_constants.h"
#include "Partons.h"
#include <iostream>
#include <string>

#include "TMath.h"
#include "TROOT.h"

/*
 *This segment of code is for calling the fortran functions related
 *to the pdf generation and initialization using the cernlibs. 
 * */

extern "C" {
    void init_pdf_(int&, int&);

    void exec_structm_(double&, double&, double&, double&, double&, double&,
                       double&, double&, double&, double&, double&);

    void exec_pkhff_(int&, int&, double&, double&, double*, double*, double*,
                     double*, double*, double*);
}

namespace HapradUtils {
    using namespace TMath;

    Double_t h3(Double_t X, Double_t q2, Double_t Z)
    {
        /* InitialiZed data */
        Double_t q0 = 1.;
        Double_t lambda = .25;
        Double_t a = -3.6544e-4;
        Double_t a1 = -2.1855;
        Double_t a2 = 3.4176;
        Double_t b1 = -1.7567;
        Double_t b2 = 1.1272;
        Double_t bb = 8.9985;
        if (q2 > q0) {
            return a * Power(X, a1) * Power((1 - X), a2) * Power(Z, b1) * Power((1 - Z), b2) * Power((Log(q2 / Power(lambda, 2))) / Power(Log(q0 / Power(lambda, 2)), 2), bb);
        } else {
            return a * Power(X, a1) * Power((1 - X), a2) * Power(Z, b1) * Power((1 - Z), b2);
        }

    }

    Double_t h4(Double_t X, Double_t q2, Double_t Z)
    {
        /* InitialiZed data */
        Double_t q0 = 1.;
        Double_t lambda = .25;
        Double_t a = .0010908;
        Double_t a1 = -3.5265e-7;
        Double_t a2 = 3.0276e-8;
        Double_t b1 = -.66787;
        Double_t b2 = 3.5868;
        Double_t bb = 6.8777;

        if (q2 > q0) {
            return a * Power(X, a1) * Power((1. - X), a2) * Power(Z, b1) * Power((1. - Z), b2) * Power(Log(q2 / Power(lambda, 2)) / Log(q0 / Power(lambda, 2)), (bb / X));
        } else {
            return a * Power(X, a1) * Power((1. - X), a2) * Power(Z, b1) * Power((1. - Z), b2);
        }
    }



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
        Double_t GTMD, SCALE, UPV, DNV, USEA, DSEA, STR, CHM, BOT, TOP, XD, GL;
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


    Double_t dfint(Int_t narg, double *arg, Int_t *nent, Double_t *ent, Double_t *table)
    {
        Double_t kd = 1;
        Double_t m = 1;
        Int_t ncomb[narg];
        Int_t ja = 0;
        Int_t jb, k;
        Int_t i = 0;
        Double_t d[narg];
        Double_t dfint;
        while (i <= narg) {
            ncomb[i] = 1;
            jb = ja + nent[i]; 
            Int_t j = ja;
            Int_t jr;
            while (j <= jb) {
                if (arg[i] <= ent[j]) {
                    if (j != ja) {
                        jr = j - 1;
                        d[i] = (ent[j] - arg[i]) / (ent[j] - ent[jr]);
                        ent[i] = j - ja;
                        kd = kd + ent[i] * m;
                        m = m * nent[i];
                        ja = jb + 1;
                        break;

                    } else {
                        j++;
                        jr = j - 1;
                        d[i] = (ent[j] - arg[i]) / (ent[j] - ent[jr]);
                        ent[i] = j - ja;
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
            ent[i] = j - ja;
            kd = kd + ent[i] * m;
            m = m * nent[i];
            ja = jb + 1;
            i++;
        }
        dfint = 0;
        while (1) {
            Double_t fac = 1.;
            Int_t iadr = kd;
            Double_t ifadr = 1.;
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
            dfint = dfint + fac * table[iadr];
            Int_t il = narg;

            while (ncomb[il] == 0) {
                il--;
                if (il == 0) return dfint;
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
        return dfint;

    }

    void ExclusiveModel(Double_t q2m, Double_t wm, Double_t csthcm, Double_t st,
                        Double_t sl, Double_t stt, Double_t stl, Double_t stlp)
    {
        const Int_t nq = 18;
        const Int_t nw = 47;
        const Int_t nt = 61;

        Double_t q2_pn[nq] = {0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1,
                              2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.2, 4.5,
                              4.8, 5.0 };

        Double_t w_pn[nw] = {1.08, 1.10, 1.12, 1.14, 1.16, 1.18, 1.20,
                             1.22, 1.24, 1.26, 1.28, 1.30, 1.32, 1.34,
                             1.36, 1.38, 1.40, 1.42, 1.44, 1.46, 1.48,
                             1.50, 1.52, 1.54, 1.56, 1.58, 1.60, 1.62,
                             1.64, 1.66, 1.68, 1.70, 1.72, 1.74, 1.76,
                             1.78, 1.80, 1.82, 1.84, 1.86, 1.88, 1.90,
                             1.92, 1.94, 1.96, 1.98, 2.00};

        Double_t th_cm_pn[nt] = {0., 3., 6., 9., 12., 15., 18., 21.,
                                 24., 27., 30., 33., 36., 39., 42.,
                                 45., 48., 51., 54., 57., 60., 63.,
                                 66., 69., 72., 75., 78., 81., 84.,
                                 87., 90., 93., 96., 99., 102., 105.,
                                 108., 111., 114., 117., 120., 123.,
                                 126., 129., 132., 135., 138., 141.,
                                 144., 147., 150., 153., 156., 159.,
                                 162., 165., 168., 171., 174., 177.,
                                 180.
                                };

        Int_t narg[3] = {nq, nw, nt};
        Double_t nc = 0;
        Double_t degrad = 57.29577952;
        Double_t a2 = 1.15;
        Double_t a3;
        Double_t a30 = -1.23;
        Double_t a31 = 0.16;
        Double_t q2cor;
        Double_t wcor;
        Double_t th_cm;
// Init
        st = 0.0;
        sl = 0.0;
        stt = 0.0;
        stl = 0.0;
        stlp = 0.0;

//new variables
        q2 = q2m;
        w = wm;
        th_cm = TMath::ACos(csthcm) * degrad;
//Check Kinematics
        if (q2 < 0.0) {
            std::cout << "Warning: Q2 < 0 in exclusive model!" << std::endl;
            std::cout << " Using Q2 = 0" << std::endl;
            q2 = 0;

        }
        if (q2 > 5) {
            q2cor = TMath::Power(5., a2) / TMath::Power(q2, a2);
            q2 = 5.;

        } else q2cor = 1.;

        if (w < 1.07784) return;

        if (w > 2) {
            if (th_cm < 50.0) a3 = a30 + a31 * 50.;
            if (th_cm > 100.0) a3 = a30 + a31 * 100.;

            a3 = a30 + a31 * th_cm;
            wcor = TMath::Power(2., a3) / TMath::Power(w, a3);
            w = 2.;

        } else wcor = 1.;

        if (TMath::Abs(csthcm) > 1) return;

        Double_t rarg[nq + nw + nt];
        const Int_t N = 100000; //solo para prueba
        double ft_cs[N], fl_cs[N],ftt_cs[N], ftl_cs[N], ftlp_cs[N];

        if (nc == 0) {
            Int_t i = 0;
            FILE *input = fopen("pi_n_maid.dat", "r");

            while (fscanf(input, "%lf   %lf   %lf   %lf   %lf\n", &ft_cs[i],&fl_cs[i], &ftt_cs[i], &ftl_cs[i], &ftlp_cs[i])) {
                i++;
            }
            fclose(input);

            for (Int_t i = 0; i < nq; i++)
                rarg[i] = q2_pn[i];

            for (Int_t j = 0; j < nw; j++)
                rarg[j + nq] = w_pn[j];

            for (Int_t k = 0; k < nt; k++)
                rarg[k + nq + nw] = th_cm_pn[k];

            nc++;

        }
        Double_t arg[3] = {q2,w,th_cm};

        st = dfint(3, arg, narg, rarg, ft_cs) * wcor * q2cor;
        sl = dfint(3, arg, narg, rarg, fl_cs) * wcor * q2cor;
        stt = dfint(3, arg, narg, rarg, ftt_cs) * wcor * q2cor;
        stl = 2. * dfint(3, arg, narg, rarg, ftl_cs) * wcor * q2cor;
        stlp = 0;

        return;
    }

}
