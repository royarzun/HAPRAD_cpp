#ifndef THAPRADUTILS_H
#define THAPRADUTILS_H

#include "TROOT.h"

namespace HapradUtils {

    Double_t fspen(const Double_t x);
    Double_t fspens(const Double_t x);

    void SemiInclusiveModel(Double_t& tldq2, Double_t& aks,
                                Double_t& tldy, Double_t& zh,
                                Double_t& tldpt2, Double_t& tldpx2,
                                Double_t& tldplh, Double_t& H1z,
                                Double_t& H2z, Double_t& H3z,
                                Double_t& H4z);
}

#endif
