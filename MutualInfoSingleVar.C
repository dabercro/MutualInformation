#include <iostream>
#include <string>
#include <TObject.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TH1.h"
#include "TMath.h"

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

void MutualInfoSingleVar(TString theFormula = "fjet1MassSDb2",TString histFileName = "fjet1MassSDb2",
                         int numBins = 80,Double_t histMin = -20,Double_t histMax = 300){

  TFile *signalFile = new TFile("signal_word.root");
  TTree *signalTree = (TTree*) signalFile->FindObjectAny("DMSTree");
  TFile *backgdFile = new TFile("background_word.root");
  TTree *backgdTree = (TTree*) backgdFile->FindObjectAny("DMSTree");

  TH1F *signalHist = new TH1F("Signal","Signal",numBins,histMin,histMax); signalHist->Sumw2();
  TH1F *backgdHist = new TH1F("Backgd","Backgd",numBins,histMin,histMax); backgdHist->Sumw2();
  TH1F *sumHist = new TH1F("Sum","Sum",numBins,histMin,histMax); sumHist->Sumw2();

  signalTree->Draw(theFormula + TString(">>Signal"),"weight*(abs(fjet1PartonId)==24)");
  backgdTree->Draw(theFormula + TString(">>Backgd"),"weight");
  *sumHist = *signalHist + *backgdHist;

  Double_t signalErr = 0.;
  Double_t signalWeights = signalHist->IntegralAndError(0,numBins,signalErr);
  fprintf(stderr,"sig int %f\n",signalWeights);
  Double_t sumErr = 0.;
  Double_t sumWeights = sumHist->IntegralAndError(0,numBins,sumErr);
  fprintf(stderr,"sum int %f\n",sumWeights);

  Double_t signalFrac = signalWeights/sumWeights;
  Double_t signalFracErr = sqrt(TMath::Power(signalErr/sumWeights,2) + TMath::Power(signalWeights*sumErr/(TMath::Power(sumWeights,2)),2));

  Double_t truthEntropy = -1*signalFrac*TMath::Log2(signalFrac) - (1-signalFrac)*TMath::Log2(1-signalFrac);
  Double_t truthErr = (TMath::Abs(1+TMath::Log2(signalFrac)) + TMath::Abs(1+TMath::Log2(1-signalFrac)))*signalFracErr;

  Double_t varErr = 0.;
  Double_t varEntropy = histEntropy(sumHist,numBins,sumWeights,sumErr,varErr);
  Double_t unionErr = 0.;
  Double_t unionEntropy = histEntropy(signalHist,numBins,sumWeights,sumErr,unionErr) +
                          histEntropy(backgdHist,numBins,sumWeights,sumErr,unionErr);

  Double_t MutualInfo = truthEntropy + varEntropy - unionEntropy;
  Double_t MutualErr = sqrt(TMath::Power(truthErr,2) + TMath::Power(varErr,2) + TMath::Power(unionErr,2));

  cout << theFormula << endl;
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
