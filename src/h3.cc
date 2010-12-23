#include "TMath.h"
#include <algorithm>
#include <cstdlib>

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
