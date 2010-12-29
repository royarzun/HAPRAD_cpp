#include "THapradUtils.h"
#include "haprad_constants.h"
#include "Partons.h"
#include <iostream>
#include <string>

#include "TMath.h"
#include "TROOT.h"


namespace HapradUtils {

    using namespace TMath;

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
            jb = ja + nent[i]; // VERIFICAR SI LOS LIMITES ESTAN CORRECTOS!!!!
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
}
