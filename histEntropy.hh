#include <iostream>
#include <string>
#include <TObject.h>

#include "TH3F.h"
#include "TH3.h"
#include "TH2F.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH1.h"

Double_t histEntropy(TH1F *aHist,int numBins,Double_t normWeight,Double_t normErr,Double_t& entropyErr){

  Double_t tempProb = 0.;
  Double_t tempErr = 0.;
  Double_t entropy = 0.;

  for(int i0 = 1;i0 <= numBins;i0++){
    tempProb = aHist->GetBinContent(i0)/normWeight;
    tempErr = sqrt(TMath::Power(aHist->GetBinError(i0)/normWeight,2) + TMath::Power(tempProb*normErr/normWeight,2));
    if(tempProb != 0){
      entropy = entropy - tempProb*TMath::Log2(tempProb);
      entropyErr = sqrt(TMath::Power(entropyErr,2) + TMath::Power(TMath::Abs(1+TMath::Log2(tempProb))*tempErr,2));
    }
  }

  return entropy;

}

Double_t histEntropy(TH2F *aHist,int numBins_x,int numBins_y,Double_t normWeight,Double_t normErr,Double_t& entropyErr){

  Double_t tempProb = 0.;
  Double_t tempErr = 0.;
  Double_t entropy = 0.;

  for(int i0 = 1;i0 < numBins_x+1;i0++){
    for(int i1 = 1;i1 < numBins_y+1;i1++){
      tempProb = aHist->GetBinContent(i0,i1)/normWeight;
      tempErr = sqrt(TMath::Power(aHist->GetBinError(i0,i1)/normWeight,2) + TMath::Power(tempProb*normErr/normWeight,2));
      if(tempProb != 0){
        entropy = entropy - tempProb*TMath::Log2(tempProb);
        entropyErr = sqrt(TMath::Power(entropyErr,2) + TMath::Power(TMath::Abs(1+TMath::Log2(tempProb))*tempErr,2));
      }
    }
  }

  return entropy;

}

Double_t histEntropy(TH3F *aHist,int numBins_x,int numBins_y,int numBins_z,Double_t normWeight,Double_t normErr,Double_t& entropyErr){

  Double_t tempProb = 0.;
  Double_t tempErr = 0.;
  Double_t entropy = 0.;

  for(int i0 = 1;i0 < numBins_x+1;i0++){
    for(int i1 = 1;i1 < numBins_y+1;i1++){
      for(int i2 = 1;i2 < numBins_z+1;i2++){
        tempProb = aHist->GetBinContent(i0,i1,i2)/normWeight;
        tempErr = sqrt(TMath::Power(aHist->GetBinError(i0,i1,i2)/normWeight,2) + TMath::Power(tempProb*normErr/normWeight,2));
        if(tempProb != 0){
          entropy = entropy - tempProb*TMath::Log2(tempProb);
          entropyErr = sqrt(TMath::Power(entropyErr,2) + TMath::Power(TMath::Abs(1+TMath::Log2(tempProb))*tempErr,2));
        }
      }
    }
  }

  return entropy;

}

