#include <iostream>
#include <string>
#include <TObject.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH3F.h"
#include "TH3.h"
#include "TH2F.h"
#include "TH2.h"
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

void MutualInfoTripleVar(TString theFormula_x = "fjet1MassSDb2",TString theFormula_y = "fjet1Tau2/fjet1Tau1",
                         TString theFormula_z = "fjet1QJetVol",TString histFileName = "MassTrimmed_Tau2Tau1",
                         int numBins_x = 80,Double_t histMin_x = -20 ,Double_t histMax_x = 300,
                         int numBins_y = 50,Double_t histMin_y = 0   ,Double_t histMax_y = 1  ,
                         int numBins_z = 60,Double_t histMin_z = -0.5,Double_t histMax_z = 2.5){

  TFile *signalFile = new TFile("signal_word.root");
  TTree *signalTree = (TTree*) signalFile->FindObjectAny("DMSTree");
  TFile *backgdFile = new TFile("background_word.root");
  TTree *backgdTree = (TTree*) backgdFile->FindObjectAny("DMSTree");

  TH3F *signalHist = new TH3F("Signal",TString("Signal;")+theFormula_x+TString(";")+theFormula_y+TString(";")+theFormula_z,
                              numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y,numBins_z,histMin_z,histMax_z);
  TH3F *backgdHist = new TH3F("Backgd",TString("Backgd;")+theFormula_x+TString(";")+theFormula_y+TString(";")+theFormula_z,
                              numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y,numBins_z,histMin_z,histMax_z);
  TH3F *sumHist = new TH3F("Sum",TString("Sum;")+theFormula_x+TString(";")+theFormula_y+TString(";")+theFormula_z,
                           numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y,numBins_z,histMin_z,histMax_z);
  signalHist->Sumw2();
  backgdHist->Sumw2();
  sumHist->Sumw2();

  signalTree->Draw(theFormula_z + TString(":") + theFormula_y + TString(":") + theFormula_x +
                   TString(">>Signal"),"weight*(abs(fjet1PartonId)==24)");
  backgdTree->Draw(theFormula_z + TString(":") + theFormula_y + TString(":") + theFormula_x +
                   TString(">>Backgd"),"weight");
  *sumHist = *signalHist + *backgdHist;

  Double_t signalErr = 0.;
  Double_t signalWeights = signalHist->IntegralAndError(1,numBins_x,1,numBins_y,1,numBins_z,signalErr);
  Double_t sumErr = 0.;
  Double_t sumWeights = sumHist->IntegralAndError(1,numBins_x,1,numBins_y,1,numBins_z,sumErr);

  Double_t signalFrac = signalWeights/sumWeights;
  Double_t signalFracErr = sqrt(TMath::Power(signalErr/sumWeights,2) + TMath::Power(signalWeights*sumErr/(TMath::Power(sumWeights,2)),2));

  Double_t truthEntropy = -1*signalFrac*TMath::Log2(signalFrac) - (1-signalFrac)*TMath::Log2(1-signalFrac);
  Double_t truthErr = (TMath::Abs(1+TMath::Log2(signalFrac)) + TMath::Abs(1+TMath::Log2(1-signalFrac)))*signalFracErr;

  Double_t varErr = 0.;
  Double_t varEntropy = histEntropy(sumHist,numBins_x,numBins_y,numBins_z,sumWeights,sumErr,varErr);
  Double_t unionErr = 0.;
  Double_t unionEntropy = histEntropy(signalHist,numBins_x,numBins_y,numBins_z,sumWeights,sumErr,unionErr) +
                          histEntropy(backgdHist,numBins_x,numBins_y,numBins_z,sumWeights,sumErr,unionErr);

  Double_t MutualInfo = truthEntropy + varEntropy - unionEntropy;
  Double_t MutualErr = sqrt(TMath::Power(truthErr,2) + TMath::Power(varErr,2) + TMath::Power(unionErr,2));

  cout << endl;
  cout << theFormula_x << " and " << theFormula_y << " and " << theFormula_z << endl;
  cout << "H(T):   " << truthEntropy << " +- " << truthErr  << endl;
  cout << "H(A):   " << varEntropy   << " +- " << varErr    << endl;
  cout << "H(T,A): " << unionEntropy << " +- " << unionErr  << endl;
  cout << "I(T;A): " << MutualInfo   << " +- " << MutualErr << endl;

  TH1F *signalHist_x = new TH1F("Signal_x","Signal",numBins_x,histMin_x,histMax_x); signalHist_x->Sumw2();
  TH1F *backgdHist_x = new TH1F("Backgd_x","Backgd",numBins_x,histMin_x,histMax_x); backgdHist_x->Sumw2();
  TH1F *signalHist_y = new TH1F("Signal_y","Signal",numBins_y,histMin_y,histMax_y); signalHist_y->Sumw2();
  TH1F *backgdHist_y = new TH1F("Backgd_y","Backgd",numBins_y,histMin_y,histMax_y); backgdHist_y->Sumw2();
  TH1F *signalHist_z = new TH1F("Signal_z","Signal",numBins_z,histMin_z,histMax_z); signalHist_z->Sumw2();
  TH1F *backgdHist_z = new TH1F("Backgd_z","Backgd",numBins_z,histMin_z,histMax_z); backgdHist_z->Sumw2();

  TH2F *signalHist_xy = new TH2F("Signal_xy","Signal",numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y); signalHist_xy->Sumw2();
  TH2F *backgdHist_xy = new TH2F("Backgd_xy","Backgd",numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y); backgdHist_xy->Sumw2();
  TH2F *signalHist_xz = new TH2F("Signal_xz","Signal",numBins_x,histMin_x,histMax_x,numBins_z,histMin_z,histMax_z); signalHist_xz->Sumw2();
  TH2F *backgdHist_xz = new TH2F("Backgd_xz","Backgd",numBins_x,histMin_x,histMax_x,numBins_z,histMin_z,histMax_z); backgdHist_xz->Sumw2();
  TH2F *signalHist_yz = new TH2F("Signal_yz","Signal",numBins_y,histMin_y,histMax_y,numBins_z,histMin_z,histMax_z); signalHist_yz->Sumw2();
  TH2F *backgdHist_yz = new TH2F("Backgd_yz","Backgd",numBins_y,histMin_y,histMax_y,numBins_z,histMin_z,histMax_z); backgdHist_yz->Sumw2();

  signalTree->Draw(theFormula_x + TString(">>Signal_x"),"weight*(abs(fjet1PartonId)==24)");
  backgdTree->Draw(theFormula_x + TString(">>Backgd_x"),"weight");
  signalTree->Draw(theFormula_y + TString(">>Signal_y"),"weight*(abs(fjet1PartonId)==24)");
  backgdTree->Draw(theFormula_y + TString(">>Backgd_y"),"weight");
  signalTree->Draw(theFormula_z + TString(">>Signal_z"),"weight*(abs(fjet1PartonId)==24)");
  backgdTree->Draw(theFormula_z + TString(">>Backgd_z"),"weight");

  signalTree->Draw(theFormula_y + TString(":") + theFormula_x + TString(">>Signal_xy"),"weight*(abs(fjet1PartonId)==24)");
  backgdTree->Draw(theFormula_y + TString(":") + theFormula_x + TString(">>Backgd_xy"),"weight");
  signalTree->Draw(theFormula_z + TString(":") + theFormula_x + TString(">>Signal_xz"),"weight*(abs(fjet1PartonId)==24)");
  backgdTree->Draw(theFormula_z + TString(":") + theFormula_x + TString(">>Backgd_xz"),"weight");
  signalTree->Draw(theFormula_z + TString(":") + theFormula_y + TString(">>Signal_yz"),"weight*(abs(fjet1PartonId)==24)");
  backgdTree->Draw(theFormula_z + TString(":") + theFormula_y + TString(">>Backgd_yz"),"weight");

  TH1F *sumHist_x = new TH1F("Sum_x","Sum",numBins_x,histMin_x,histMax_x); sumHist_x->Sumw2();
  TH1F *sumHist_y = new TH1F("Sum_y","Sum",numBins_y,histMin_y,histMax_y); sumHist_y->Sumw2();
  TH1F *sumHist_z = new TH1F("Sum_z","Sum",numBins_z,histMin_z,histMax_z); sumHist_z->Sumw2();
  TH2F *sumHist_xy = new TH2F("Sum_xy","Sum",numBins_x,histMin_x,histMax_x,numBins_y,histMin_y,histMax_y); sumHist_xy->Sumw2();
  TH2F *sumHist_xz = new TH2F("Sum_xz","Sum",numBins_x,histMin_x,histMax_x,numBins_z,histMin_z,histMax_z); sumHist_xz->Sumw2();
  TH2F *sumHist_yz = new TH2F("Sum_yz","Sum",numBins_y,histMin_y,histMax_y,numBins_z,histMin_z,histMax_z); sumHist_yz->Sumw2();

  *sumHist_x = *signalHist_x + *backgdHist_x;
  *sumHist_y = *signalHist_y + *backgdHist_y;
  *sumHist_z = *signalHist_z + *backgdHist_z;
  *sumHist_xy = *signalHist_xy + *backgdHist_xy;
  *sumHist_xz = *signalHist_xz + *backgdHist_xz;
  *sumHist_yz = *signalHist_yz + *backgdHist_yz;

  Double_t sumErr_x = 0.;
  Double_t sum_x = sumHist_x->IntegralAndError(1,numBins_x,sumErr_x);
  Double_t sumErr_y = 0.;
  Double_t sum_y = sumHist_y->IntegralAndError(1,numBins_y,sumErr_y);
  Double_t sumErr_z = 0.;
  Double_t sum_z = sumHist_z->IntegralAndError(1,numBins_z,sumErr_z);
  Double_t sumErr_xy = 0.;
  Double_t sum_xy = sumHist_xy->IntegralAndError(1,numBins_x,1,numBins_y,sumErr_xy);
  Double_t sumErr_xz = 0.;
  Double_t sum_xz = sumHist_xz->IntegralAndError(1,numBins_x,1,numBins_z,sumErr_xz);
  Double_t sumErr_yz = 0.;
  Double_t sum_yz = sumHist_yz->IntegralAndError(1,numBins_y,1,numBins_z,sumErr_yz);

  Double_t xErr = 0.;
  Double_t xEntropy = histEntropy(sumHist_x,numBins_x,sum_x,sumErr_x,xErr);
  Double_t yErr = 0.;
  Double_t yEntropy = histEntropy(sumHist_y,numBins_y,sum_y,sumErr_y,yErr);
  Double_t zErr = 0.;
  Double_t zEntropy = histEntropy(sumHist_z,numBins_z,sum_z,sumErr_z,zErr);
  Double_t xyErr = 0.;
  Double_t xyEntropy = histEntropy(sumHist_xy,numBins_x,numBins_y,sum_xy,sumErr_xy,xyErr);
  Double_t xzErr = 0.;
  Double_t xzEntropy = histEntropy(sumHist_xz,numBins_x,numBins_z,sum_xz,sumErr_xz,xzErr);
  Double_t yzErr = 0.;
  Double_t yzEntropy = histEntropy(sumHist_yz,numBins_y,numBins_z,sum_yz,sumErr_yz,yzErr);

  Double_t Mutual_xy = xEntropy + yEntropy - xyEntropy;
  Double_t Mutual_xz = xEntropy + zEntropy - xzEntropy;
  Double_t Mutual_yz = yEntropy + zEntropy - yzEntropy;
  Double_t MutualErr_xy = sqrt(TMath::Power(xErr,2) + TMath::Power(yErr,2) + TMath::Power(xyErr,2));
  Double_t MutualErr_xz = sqrt(TMath::Power(xErr,2) + TMath::Power(zErr,2) + TMath::Power(xzErr,2));
  Double_t MutualErr_yz = sqrt(TMath::Power(yErr,2) + TMath::Power(zErr,2) + TMath::Power(yzErr,2));

  cout << endl;
  cout << "Correlation Information" << endl;
  cout << "A:      " << theFormula_x << endl;
  cout << "B:      " << theFormula_y << endl;
  cout << "C:      " << theFormula_z << endl;
  cout << "H(A):   " << xEntropy     << " +- " << xErr  << endl;
  cout << "H(B):   " << yEntropy     << " +- " << yErr  << endl;
  cout << "H(C):   " << zEntropy     << " +- " << zErr  << endl;
  cout << "H(A,B): " << xyEntropy    << " +- " << xyErr << endl;
  cout << "H(A,C): " << xzEntropy    << " +- " << xzErr << endl;
  cout << "H(B,C): " << yzEntropy    << " +- " << yzErr << endl;
  cout << "I(A,B): " << Mutual_xy    << " +- " << MutualErr_xy << endl;
  cout << "I(A,C): " << Mutual_xz    << " +- " << MutualErr_xz << endl;
  cout << "I(B,C): " << Mutual_yz    << " +- " << MutualErr_yz << endl;
  cout << endl;

  TFile *outFile = new TFile(TString("scratch/") + histFileName + TString(".root"),"RECREATE");

  signalHist->Write();
  backgdHist->Write();
  sumHist->Write();

  signalHist_x->Write();
  backgdHist_x->Write();
  sumHist_x->Write();
  signalHist_y->Write();
  backgdHist_y->Write();
  sumHist_y->Write();
  signalHist_z->Write();
  backgdHist_z->Write();
  sumHist_z->Write();
  signalHist_xy->Write();
  backgdHist_xy->Write();
  sumHist_xy->Write();
  signalHist_xz->Write();
  backgdHist_xz->Write();
  sumHist_xz->Write();
  signalHist_yz->Write();
  backgdHist_yz->Write();
  sumHist_yz->Write();

  outFile->Close();

}
