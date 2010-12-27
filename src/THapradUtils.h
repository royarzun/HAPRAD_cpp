#ifndef THAPRADUTILS_H
#define THAPRADUTILS_H

#include "TROOT.h"

namespace HapradUtils {

    Double_t fspen(const Double_t x);
    Double_t fspens(const Double_t x);
    void SemiInclusiveModel(Double_t q2, Double_t X,
                            Double_t Y, Double_t Z,
                            Double_t pt2, Double_t mx2,
                            Double_t pl, Double_t& H1z,
                            Double_t& H2z, Double_t& H3z,
                            Double_t& H4z);

    Double_t dfint(Int_t narg, double *arg, Int_t *nent, Double_t *ent, Double_t *table);
    void ExclusiveModel(Double_t q2m, Double_t wm, Double_t csthcm, Double_t st,
                        Double_t sl, Double_t stt, Double_t stl, Double_t stlp);

}

#endif
