#include "TMath.h"

using namespace TMath;

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
		return a*Power(X,a1)*Power((1.-X),a2)*Power(Z,b1)*Power((1.-Z),b2)*Power(Log(q2/Power(lambda,2))/Log(q0/Power(lambda,2)),(bb/X));

    } else {
		return a*Power(X,a1)*Power((1.-X),a2)*Power(Z,b1)*Power((1.-Z),b2);
    }
} 

