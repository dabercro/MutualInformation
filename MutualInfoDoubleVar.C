#include <iostream>
#include <string>
#include <TObject.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH2F.h"
#include "TH2.h"
#include "TH1.h"
#include "TMath.h"

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

void MutualInfoDoubleVar(TString theFormula_x = "fjet1MassSDb2",TString theFormula_y = "fjet1Tau2/fjet1Tau1",
                         TString histFileName = "MassTrimmed_Tau2Tau1",
                         int numBins_x = 80,Double_t histMin_x = -20,Double_t histMax_x = 300,
                         int numBins_y = 50,Double_t histMin_y = 0  ,Double_t histMax_y = 1  ){

  TFile *signalFile = new TFile("signal_word.root");
  TTree *signalTree = (TTree*) signalFile->FindObjectAny("DMSTree");
  TFile *backgdFile = new TFile("background_word.root");
  TTree *backgdTree = (TTree*) backgdFile->FindObjectAny("DMSTree");

  TH2F *signalHist = new TH2F("Signal",TString("Signal;")+theFormula_x+TString(";")+theFormula_y,
                              numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y);
  TH2F *backgdHist = new TH2F("Backgd",TString("Backgd;")+theFormula_x+TString(";")+theFormula_y,
                              numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y);
  TH2F *sumHist = new TH2F("Sum",TString("Sum;")+theFormula_x+TString(";")+theFormula_y,
                           numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y);
  signalHist->Sumw2();
  backgdHist->Sumw2();
  sumHist->Sumw2();

  signalTree->Draw(theFormula_y + TString(":") + theFormula_x + TString(">>Signal"),"weight*(abs(fjet1PartonId)==24)");
  backgdTree->Draw(theFormula_y + TString(":") + theFormula_x + TString(">>Backgd"),"weight");
  *sumHist = *signalHist + *backgdHist;

  Double_t signalErr = 0.;
  Double_t signalWeights = signalHist->IntegralAndError(1,numBins_x,1,numBins_y,signalErr);
  Double_t sumErr = 0.;
  Double_t sumWeights = sumHist->IntegralAndError(1,numBins_x,1,numBins_y,sumErr);

  Double_t signalFrac = signalWeights/sumWeights;
  Double_t signalFracErr = sqrt(TMath::Power(signalErr/sumWeights,2) + TMath::Power(signalWeights*sumErr/(TMath::Power(sumWeights,2)),2));

  Double_t truthEntropy = -1*signalFrac*TMath::Log2(signalFrac) - (1-signalFrac)*TMath::Log2(1-signalFrac);
  Double_t truthErr = (TMath::Abs(1+TMath::Log2(signalFrac)) + TMath::Abs(1+TMath::Log2(1-signalFrac)))*signalFracErr;

  Double_t varErr = 0.;
  Double_t varEntropy = histEntropy(sumHist,numBins_x,numBins_y,sumWeights,sumErr,varErr);
  Double_t unionErr = 0.;
  Double_t unionEntropy = histEntropy(signalHist,numBins_x,numBins_y,sumWeights,sumErr,unionErr) +
                          histEntropy(backgdHist,numBins_x,numBins_y,sumWeights,sumErr,unionErr);

  Double_t MutualInfo = truthEntropy + varEntropy - unionEntropy;
  Double_t MutualErr = sqrt(TMath::Power(truthErr,2) + TMath::Power(varErr,2) + TMath::Power(unionErr,2));

  cout << theFormula_x << " and " << theFormula_y << endl;
  cout << "H(T):   " << truthEntropy << " +- " << truthErr  << endl;
  cout << "H(A):   " << varEntropy   << " +- " << varErr    << endl;
  cout << "H(T,A): " << unionEntropy << " +- " << unionErr  << endl;
  cout << "I(T;A): " << MutualInfo   << " +- " << MutualErr << endl;

  TFile *outFile = new TFile(TString("scratch/") + histFileName + TString(".root"),"RECREATE");
  signalHist->Write();
  backgdHist->Write();
  sumHist->Write();
  outFile->Close();

}
