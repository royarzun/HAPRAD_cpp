#include "TMath.h"

const Double_t kBarn = 0.389379 * TMath::Power(10,6);
const Double_t kAlpha = 0.729735 * TMath::Power(10,-2);
const Double_t kPi = TMath::Pi();

const Double_t kMassProton = 0.938272;
const Double_t kMassNeutron = 0.93956536;
const Double_t kMassElectron = 0.000511;
const Double_t kMassMuon = 0.10565;
const Double_t kMassPion = 0.1395675;

const Double_t kMassDetectedHadron   = 0.1395675;
const Double_t kMassUndetectedHadron = 0.93956536;

const Double_t kMassC2 = TMath::Power(kMassProton + kMassDetectedHadron,2);

const Double_t kRadianDeg = 57.2957795131;
const Double_t kEpsMachine = 1.1 * TMath::Power(10,-15);
