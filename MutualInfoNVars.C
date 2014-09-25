#include <iostream>
#include <string>
#include <TObject.h>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TH1.h"
#include "TMath.h"

Double_t histEntropy(THnF *aHist,int dim, int* numBins,Double_t normWeight,Double_t normErr,Double_t& entropyErr){

  Double_t tempProb = 0.;
  Double_t tempErr = 0.;
  Double_t entropy = 0.;
  int *indices=new int[dim];
  for(int i=0;i<dim;i++) {
    indices[i]=1; //the first bin
  }
  while (1) {
    tempProb = aHist->GetBinContent(indices)/normWeight;
    tempErr = sqrt(TMath::Power(aHist->GetBinError(indices)/normweight,2) + TMath::Power(tempProb*normErr/normWeight,2));
    if (tempProb) {
      entropy = entropy - tempProb*TMath::Log2(tempProb);
      entropyErr = sqrt(TMath::Power(entropyErr,2) + TMath::Power(TMath::Abs(1+TMath::Log2(tempProb))*tempErr,2));
    }
    int j;
    for(j=0;j<dim;j++) {
      indices[j]++;
      if (indices[j] < numBins[j]) break; //this is a valid index
      else indices[j]=0; //done with valid indices, loop around
    }
    if(j==dim) break;
  }

  return entropy;

}

void MutualInfoNVars(TString theFormula = "fjet1MassSDb2",TString histFileName = "fjet1MassSDb2", Int_t dim=1,Int_t numBin[] = 0,Float_t histMin[] = 0,Float_t histMax[] = 0){
    if (dim==1) { // because vim syntax highlighting is stupid don't ask
        if (numBin==0) {
            numBin = malloc(1*sizeof(Int_t));
            numBin[0]=80;
        }
        if (histMin==0) {
            histMin = malloc(1*sizeof(Float_t));
            histMin[0]=-20;
        }
        if (histMax==0) {
            histMax = malloc(1*sizeof(Float_t));
            histMax[0]=300;
        }
    }

    Int_t expMem=1;
    for(int i=0;i<dim;i++)  expMem*=numBin[i];
    expMem*=6*sizeof(Double_t);
    //don't use so much memory you break hadoop
    fprintf(stderr,"WARNING: You are about to use %i B of memory\n",expMem);
    fprintf(stderr,"Pausing for %.2f seconds so you can re-evaluate and hit Ctrl-C",float(expMem)/(1024*1024*512));
    usleep(expMem*100/(1024*1024*512));
    
  TFile *signalFile = new TFile("signal_word.root");
  TTree *signalTree = (TTree*) signalFile->FindObjectAny("DMSTree");
  TFile *backgdFile = new TFile("background_word.root");
  TTree *backgdTree = (TTree*) backgdFile->FindObjectAny("DMSTree");

  THnD *signalHist = new THnD("Signal","Signal",dim,numBin,histMin,histMax);  signalHist->Sumw2();
  THnD *backgdHist = new THnD("Backgd","Backgd",dim,numBin,histMin,histMax);  backgdHist->Sumw2();
  THnD *sumHist ; //= new THnF("Sum","Sum",numBins,histMin,histMax); sumHist->Sumw2();
  
  signalTree->Draw(theFormula + TString(">>Signal"),"weight*(abs(fjet1PartonId)==24)");
  backgdTree->Draw(theFormula + TString(">>Backgd"),"weight");
  sumHist = backgdHist->Clone("Total");
  sumHist->Add(signalHist); 
  Double_t signalErr = 0.;
  Double_t signalIntegral = signalHist->IntegralAndError(1,numBins,signalErr);
  Double_t sumErr = 0.;
  Double_t sumIntegral = sumHist->IntegralAndError(1,numBins,sumErr);

  Double_t signalFrac = signalIntegral/sumIntegral;
  Double_t signalFracErr = sqrt(TMath::Power(signalErr/sumIntegral,2) + TMath::Power(signalIntegral*sumErr/(TMath::Power(sumIntegral,2)),2));

  Double_t truthEntropy = -1*signalFrac*TMath::Log2(signalFrac) - (1-signalFrac)*TMath::Log2(1-signalFrac);
  Double_t truthErr = (TMath::Abs(1+TMath::Log2(signalFrac)) + TMath::Abs(1+TMath::Log2(1-signalFrac)))*signalFracErr;

  Double_t varErr = 0.;
  Double_t varEntropy = histEntropy(sumHist,numBins,sumIntegral,sumErr,varErr);
  Double_t unionErr = 0.;
  Double_t unionEntropy = histEntropy(signalHist,numBins,sumIntegral,sumErr,unionErr) +
                          histEntropy(backgdHist,numBins,sumIntegral,sumErr,unionErr);

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
