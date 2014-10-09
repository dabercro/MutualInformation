#include <iostream>
#include <string>
#include <TObject.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH2F.h"
#include "TH2.h"
#include "TMath.h"

#include "histEntropy.hh"

/////////////////////////////////////////////////////////////////////
// Note: Bias reduction hasn't been demonstrated with this yet.    //
//       At least, the uncertainty seems larger than current bias. //
//       I will try to get the bias to show itself though.         //
//       Perhaps the higher-dimensional cases will do this.        //
//       Also, I could just be doing something wrong...            //
/////////////////////////////////////////////////////////////////////

void MutualInfoDoubleVar_reducedBias(TString theFormula_x = "fjet1MassTrimmed",TString theFormula_y = "fjet1PullAngle",
                         TString histFileName = "test",
                         int numBins_x = 210,Double_t histMin_x =  -20,Double_t histMax_x = 400,
                         int numBins_y = 140,Double_t histMin_y = -3.5,Double_t histMax_y = 3.5){

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
  Double_t backgdErr = 0.;
  Double_t backgdWeights = backgdHist->IntegralAndError(1,numBins_x,1,numBins_y,backgdErr);
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

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //                            Attempts to reduce bias with A.20                                        //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  Double_t signalEntropyErr = 0.;
  Double_t signalEntropy = histEntropy(signalHist,numBins_x,numBins_y,signalWeights,signalErr,signalEntropyErr);
  Double_t backgdEntropyErr = 0.;
  Double_t backgdEntropy = histEntropy(backgdHist,numBins_x,numBins_y,backgdWeights,backgdErr,backgdEntropyErr);

  TH2F *signalHist_0 = new TH2F("Signal_0",TString("Signal_0;")+theFormula_x+TString(";")+theFormula_y,
                              numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y);
  TH2F *backgdHist_0 = new TH2F("Backgd_0",TString("Backgd_0;")+theFormula_x+TString(";")+theFormula_y,
                              numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y);
  TH2F *sumHist_0 = new TH2F("Sum_0",TString("Sum_0;")+theFormula_x+TString(";")+theFormula_y,
                           numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y);

  signalHist_0->Sumw2();
  backgdHist_0->Sumw2();
  sumHist_0->Sumw2();

  // This event % 2 term splits the sample in half
  signalTree->Draw(theFormula_y + TString(":") + theFormula_x + TString(">>Signal_0"),"weight*(abs(fjet1PartonId)==24 && event % 2 == 0)");
  backgdTree->Draw(theFormula_y + TString(":") + theFormula_x + TString(">>Backgd_0"),"weight*(event % 2 == 0)"); 
  *sumHist_0 = *signalHist_0 + *backgdHist_0;

  Double_t sumErr_0 = 0.;
  Double_t sumWeights_0 = sumHist_0->IntegralAndError(1,numBins_x,1,numBins_y,signalErr);

  Double_t mixedErr_0 = 0.;
  Double_t mixedEntropy_0 = histEntropy(sumHist_0,numBins_x,numBins_y,sumWeights_0,sumErr_0,mixedErr_0);

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

  TH2F *signalHist_1 = new TH2F("Signal_1",TString("Signal_1;")+theFormula_x+TString(";")+theFormula_y,
                              numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y);
  TH2F *backgdHist_1 = new TH2F("Backgd_1",TString("Backgd_1;")+theFormula_x+TString(";")+theFormula_y,
                              numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y);
  TH2F *sumHist_1 = new TH2F("Sum_1",TString("Sum_1;")+theFormula_x+TString(";")+theFormula_y,
                           numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y);

  signalHist_1->Sumw2();
  backgdHist_1->Sumw2();
  sumHist_1->Sumw2();

  // This event % 2 term splits the sample in half
  signalTree->Draw(theFormula_y + TString(":") + theFormula_x + TString(">>Signal_1"),"weight*(abs(fjet1PartonId)==24 && event % 2 == 1)");
  backgdTree->Draw(theFormula_y + TString(":") + theFormula_x + TString(">>Backgd_1"),"weight*(event % 2 == 1)"); 
  *sumHist_1 = *signalHist_1 + *backgdHist_1;

  Double_t sumErr_1 = 0.;
  Double_t sumWeights_1 = sumHist_1->IntegralAndError(1,numBins_x,1,numBins_y,signalErr);

  Double_t mixedErr_1 = 0.;
  Double_t mixedEntropy_1 = histEntropy(sumHist_1,numBins_x,numBins_y,sumWeights_1,sumErr_1,mixedErr_1);

  Double_t MutualInfo_1 = mixedEntropy_1 - signalFrac*signalEntropy - (1-signalFrac)*backgdEntropy;
  Double_t MutualErr_1  = sqrt(TMath::Power(mixedErr_1,2) +
                               TMath::Power(signalFracErr*(signalEntropy-backgdEntropy),2) +
                               TMath::Power(signalFrac*signalEntropyErr,2) +
                               TMath::Power((1-signalFrac)*backgdEntropyErr,2));

  cout << endl;
  cout << "Reduced Bias Results (0)" << endl;
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
