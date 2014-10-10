#include <iostream>
#include <string>
#include <TObject.h>

#include "THnSparse.h"

Float_t IntAndErr(THnSparseF*aHist, int dim, int *numBins,Float_t &integralErr) {
//  Double_t integral=aHist->ComputeIntegral(); // I don't know why this always returns 1
//  GetBinContent/Error works, so compute manually
  Float_t integral=0.;
  integralErr = 0;
  Int_t nFilledBins = aHist->GetNbins();
  for (int i=0; i < nFilledBins; i++) {
    integral+=aHist->GetBinContent(i);
    integralErr+=TMath::Power(aHist->GetBinError(i),2.);
  }
 integralErr=sqrt(integralErr);
 return integral;   
}

