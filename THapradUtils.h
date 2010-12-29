#ifndef THAPRADUTILS_H
#define THAPRADUTILS_H

#include "TROOT.h"

namespace HapradUtils {

    Double_t fspen(const Double_t x);
    Double_t fspens(const Double_t x);
    Double_t dfint(Int_t narg, double *arg, Int_t *nent, Double_t *ent, Double_t *table);
    void ExclusiveModel(Double_t q2m, Double_t wm, Double_t csthcm, Double_t st,
                        Double_t sl, Double_t stt, Double_t stl, Double_t stlp);

}

#endif
