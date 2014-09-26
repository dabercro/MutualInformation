#include <iostream>
#include <string>
#include <TObject.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TH1.h"
#include "TMath.h"

/////////////////////////////////////////////////////////////////////
// Note: Bias reduction hasn't been demonstrated with this yet.    //
//       At least, the uncertainty seems larger than current bias. //
//       I will try to get the bias to show itself though.         //
//       Perhaps the higher-dimensional cases will do this.        //
//       Also, I could just be doing something wrong...            //
/////////////////////////////////////////////////////////////////////

Double_t histEntropy(TH1F *aHist,int numBins,Double_t normWeight,Double_t normErr,Double_t& entropyErr){

  Double_t tempProb = 0.;
  Double_t tempErr = 0.;
  Double_t entropy = 0.;

  for(int i0 = 1;i0 < numBins+1;i0++){
    tempProb = aHist->GetBinContent(i0)/normWeight;
    tempErr = sqrt(TMath::Power(aHist->GetBinError(i0)/normWeight,2) + TMath::Power(tempProb*normErr/normWeight,2));
    if(tempProb != 0){
      entropy = entropy - tempProb*TMath::Log2(tempProb);
      entropyErr = sqrt(TMath::Power(entropyErr,2) + TMath::Power(TMath::Abs(1+TMath::Log2(tempProb))*tempErr,2));
    }
  }

  return entropy;

}

void MutualInfoSingleVar_reducedBias(TString theFormula = "fjet1MassSDb2",TString histFileName = "fjet1MassSDb2",
                                     int numBins = 320,Double_t histMin = -20,Double_t histMax = 300){

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
  Double_t signalWeights = signalHist->IntegralAndError(1,numBins,signalErr);
  Double_t backgdErr = 0.;
  Double_t backgdWeights = backgdHist->IntegralAndError(1,numBins,backgdErr);
  Double_t sumErr = 0.;
  Double_t sumWeights = sumHist->IntegralAndError(1,numBins,sumErr);

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

  cout << endl;
  cout << theFormula << endl;
  cout << endl;
  cout << "Old Way" << endl;
  cout << "H(T):   " << truthEntropy << " +- " << truthErr  << endl;
  cout << "H(A):   " << varEntropy   << " +- " << varErr    << endl;
  cout << "H(T,A): " << unionEntropy << " +- " << unionErr  << endl;
  cout << "I(T;A): " << MutualInfo   << " +- " << MutualErr << endl;

  TFile *outFile = new TFile(TString("scratch/") + histFileName + TString(".root"),"RECREATE");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                            Attempts to reduce bias with A.20                                        //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  Double_t signalEntropyErr = 0.;
  Double_t signalEntropy = histEntropy(signalHist,numBins,signalWeights,signalErr,signalEntropyErr);
  Double_t backgdEntropyErr = 0.;
  Double_t backgdEntropy = histEntropy(backgdHist,numBins,backgdWeights,backgdErr,backgdEntropyErr);

  TH1F *signalHist_0 = new TH1F("Signal_0","Signal_0",numBins,histMin,histMax); signalHist_0->Sumw2();
  TH1F *backgdHist_0 = new TH1F("Backgd_0","Backgd_0",numBins,histMin,histMax); backgdHist_0->Sumw2();
  TH1F *sumHist_0 = new TH1F("Sum_0","Sum_0",numBins,histMin,histMax); sumHist_0->Sumw2();

  signalTree->Draw(theFormula + TString(">>Signal_0"),"weight*(abs(fjet1PartonId)==24 && event % 2 == 0)");
  backgdTree->Draw(theFormula + TString(">>Backgd_0"),"weight*(event % 2 == 0)"); // This event % 2 term splits the sample in half
  *sumHist_0 = *signalHist_0 + *backgdHist_0;

  Double_t sumErr_0 = 0.;
  Double_t sumWeights_0 = sumHist_0->IntegralAndError(1,numBins,sumErr);

  Double_t mixedErr_0 = 0.;
  Double_t mixedEntropy_0 = histEntropy(sumHist_0,numBins,sumWeights_0,sumErr_0,mixedErr_0);

  Double_t MutualInfo_0 = mixedEntropy_0 - signalFrac*signalEntropy - (1-signalFrac)*backgdEntropy;
  Double_t MutualErr_0  = sqrt(TMath::Power(mixedErr_0,2) +
                               TMath::Power(signalFracErr*(signalEntropy-backgdEntropy),2) +
                               TMath::Power(signalFrac*signalEntropyErr,2) +
                               TMath::Power((1-signalFrac)*backgdEntropyErr,2));

  cout << endl;
  cout << "Reduced Bias Results (0)" << endl;
  cout << "f:         " << signalFrac     << " +- " << signalFracErr    << endl;
  cout << "H_sig+bck: " << mixedEntropy_0 << " +- " << mixedErr_0       << endl;
  cout << "H_sig:     " << signalEntropy  << " +- " << signalEntropyErr << endl;
  cout << "H_bck:     " << backgdEntropy  << " +- " << backgdEntropyErr << endl;
  cout << "I(T;A):    " << MutualInfo_0   << " +- " << MutualErr_0      << endl;

  TH1F *signalHist_1 = new TH1F("Signal_1","Signal_1",numBins,histMin,histMax); signalHist_1->Sumw2();
  TH1F *backgdHist_1 = new TH1F("Backgd_1","Backgd_1",numBins,histMin,histMax); backgdHist_1->Sumw2();
  TH1F *sumHist_1 = new TH1F("Sum_1","Sum_1",numBins,histMin,histMax); sumHist_1->Sumw2();

  signalTree->Draw(theFormula + TString(">>Signal_1"),"weight*(abs(fjet1PartonId)==24 && event % 2 == 1)");
  backgdTree->Draw(theFormula + TString(">>Backgd_1"),"weight*(event % 2 == 1)"); // This event % 2 term splits the sample in half
  *sumHist_1 = *signalHist_1 + *backgdHist_1;

  Double_t sumErr_1 = 0.;
  Double_t sumWeights_1 = sumHist_1->IntegralAndError(1,numBins,sumErr);

  Double_t mixedErr_1 = 0.;
  Double_t mixedEntropy_1 = histEntropy(sumHist_1,numBins,sumWeights_1,sumErr_1,mixedErr_1);

  Double_t MutualInfo_1 = mixedEntropy_1 - signalFrac*signalEntropy - (1-signalFrac)*backgdEntropy;
  Double_t MutualErr_1  = sqrt(TMath::Power(mixedErr_1,2) +
                               TMath::Power(signalFracErr*(signalEntropy-backgdEntropy),2) +
                               TMath::Power(signalFrac*signalEntropyErr,2) +
                               TMath::Power((1-signalFrac)*backgdEntropyErr,2));

  cout << endl;
  cout << "Reduced Bias Results (1)" << endl;
  cout << "f:         " << signalFrac     << " +- " << signalFracErr    << endl;
  cout << "H_sig+bck: " << mixedEntropy_1 << " +- " << mixedErr_1       << endl;
  cout << "H_sig:     " << signalEntropy  << " +- " << signalEntropyErr << endl;
  cout << "H_bck:     " << backgdEntropy  << " +- " << backgdEntropyErr << endl;
  cout << "I(T;A):    " << MutualInfo_1   << " +- " << MutualErr_1      << endl;

  signalHist->Write();
  backgdHist->Write();
  sumHist->Write();
  signalHist_0->Write();
  backgdHist_0->Write();
  sumHist_0->Write();
  signalHist_1->Write();
  backgdHist_1->Write();
  sumHist_1->Write();

  outFile->Close();

}
