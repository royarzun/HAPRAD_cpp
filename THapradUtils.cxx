#include "THapradUtils.h"
#include "TMath.h"


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
            logprod = TMath::Log(x) * TMath::Log(1 - x + TMath::Power(10,-10));
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
        Double_t tch = TMath::Power(10,-16);
        Double_t b;

        do {
            an = an + 1;
            a = a * x;
            b = a / TMath::Power(an,2);
            f = f + b;
        } while (b - tch > 0);

        return f;
    }



    void SemiInclusiveModel(Double_t& tldq2, Double_t& aks,
                                Double_t& tldy, Double_t& zh,
                                Double_t& tldpt2, Double_t& tldpx2,
                                Double_t& tldplh, Double_t& H1z,
                                Double_t& H2z, Double_t& H3z,
                                Double_t& H4z)
    {

    }
}
