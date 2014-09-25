#include <iostream>
#include <string>
#include <TObject.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TH1.h"
#include "TMath.h"
#include "THn.h"

float histEntropy(THnF *aHist,int dim, int* numBin,float normWeight){

  float tempProb = 0.;
  float entropy = 0.;
  Int_t *
  for(int d=0; d<dim; d++) {
    for(int i0 = 0;i0 < numBin[d];i0++){
      tempProb = aHist->GetBinContent(i0)/normWeight;
      if(tempProb != 0) entropy = entropy - tempProb*TMath::Log2(tempProb);
    }
  }

  return entropy;

}

void MutualInfoNVars(TString theFormula = "fjet1MassSDb2",int dim=1,int numBin = 80,float histMin = -20,float histMax = 300){

  TFile *signalFile = new TFile("signal_word.root");
  TTree *signalTree = (TTree*) signalFile->FindObjectAny("DMSTree");
  TFile *backgdFile = new TFile("background_word.root");
  TTree *backgdTree = (TTree*) backgdFile->FindObjectAny("DMSTree");

  TH1F *signalHist = new TH1F("Signal","Signal",numBin,histMin,histMax);
  signalHist->Sumw2();
  TH1F *backgdHist = new TH1F("Backgd","Backgd",numBin,histMin,histMax);
  backgdHist->Sumw2();
  TH1F *sumHist = new TH1F("Sum","Sum",numBin,histMin,histMax);
  sumHist->Sumw2();

  // TH1F *signalHalfHist = new TH1F("SignalHalf","SignalHalf",numBin/2+1,histMin,histMax);
  // signalHalfHist->Sumw2();
  // TH1F *backgdHalfHist = new TH1F("BackgdHalf","BackgdHalf",numBin/2+1,histMin,histMax);
  // backgdHalfHist->Sumw2();

  signalTree->Draw(theFormula + TString(">>Signal"),"weight*(abs(fjet1PartonId)==24)");
  backgdTree->Draw(theFormula + TString(">>Backgd"),"weight");
  *sumHist = *signalHist + *backgdHist;
  // signalTree->Draw(theFormula + TString(">>SignalHalf"),"weight*(abs(fjet1PartonId)==24)");
  // backgdTree->Draw(theFormula + TString(">>BackgdHalf"),"weight");

  float signalWeights = signalHist->GetSumOfWeights();
  float backgdWeights = backgdHist->GetSumOfWeights();
  float sumWeights = sumHist->GetSumOfWeights();
  float signalFrac = signalWeights/(signalWeights + backgdWeights);

  float truthEntropy = -1*signalFrac*TMath::Log2(signalFrac) - (1-signalFrac)*TMath::Log2(1-signalFrac);

  float varEntropy = histEntropy(sumHist,numBin,sumWeights);
  float unionEntropy = histEntropy(signalHist,numBin,sumWeights) + histEntropy(backgdHist,numBin,sumWeights);
  // float unionEntropy = histEntropy(signalHalfHist,numBin/2+1,sumWeights) + histEntropy(backgdHalfHist,numBin/2+1,sumWeights);

  float MutualInfo = truthEntropy + varEntropy - unionEntropy;

  cout << "H(T):   " << truthEntropy << endl;
  cout << "H(A):   " << varEntropy << endl;
  cout << "H(T,A): " << unionEntropy << endl;
  cout << "I(T;A): " << MutualInfo << endl;

  TFile *outFile = new TFile(TString("scratch/") + theFormula + TString(".root"),"RECREATE");
  signalHist->Write();
  backgdHist->Write();
  sumHist->Write();
  outFile->Close();

}
