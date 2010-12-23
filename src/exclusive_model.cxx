#include <iostream>
#include <fstream.h>
#include <string>
#include "TMath.h"
#include "TRoot.h"

#define N 100000;

void exclusive_model(Double_t q2m, Double_t wm, Double_t csthcm, Double_t st, Double_t sl, Double_t stt, Double_t stl, Double_t stlp)
{
    Double_t nq = 18;
    Double_t nw = 47;
    Double_t nt = 61;

    Double_t q2_pn[nq] = {0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8,
                          2.1, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9,
                          4.2, 4.5, 4.8, 5.0
                         };

    Double_t w_pn[nw] = {1.08, 1.10, 1.12, 1.14, 1.16, 1.18,
                         1.20, 1.22, 1.24, 1.26, 1.28, 1.30,
                         1.32, 1.34, 1.36, 1.38, 1.40, 1.42,
                         1.44, 1.46, 1.48, 1.50, 1.52, 1.54,
                         1.56, 1.58, 1.60, 1.62, 1.64, 1.66,
                         1.68, 1.70, 1.72, 1.74, 1.76, 1.78,
                         1.80, 1.82, 1.84, 1.86, 1.88, 1.90,
                         1.92, 1.94, 1.96, 1.98, 2.00
                        };
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
    Double_t a30 = -1.23;
    Double_t a31 = 0.16;
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
    Double_t ft_cs[N], fl_cs[N], ftt_cs[N], ftl_cs[N], ftlp_cs[N];

    if (nc == 0) {
        Int_t i = 0;
        FILE *input = fopen("pi_n_maid.dat", "r");

        while (fscanf(input, "%f   %f   %f   %f   %f", &ft_cs[i], &fl_cs[i], &ftt_cs[i], &ftl_cs, &ftlp_cs[i])) {
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

    arg[0] = q2;
    arg[1] = w;
    arg[2] = th_cm;

    st = dfint(3, arg, narg, rarg, ft_cs) * wcor * q2cor;
    sl = dfint(3, arg, narg, rarg, fl_cs) * wcor * q2cor;
    stt = dfint(3, arg, narg, rarg, ftt_cs) * wcor * q2cor;
    stl = 2. * dfint(3, arg, narg, rarg, ftl_cs) * wcor * q2cor;
    stlp = 0;

	return;
}
