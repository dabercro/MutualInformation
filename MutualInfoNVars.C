#include <iostream>
#include <fstream>
#include <string>
#include <TObject.h>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TH1.h"
#include "TMath.h"
#include "THn.h"

Float_t histEntropy(THnSparseF *aHist,int dim, int* numBins,Float_t normWeight,Float_t normErr,Float_t& entropyErr){

  Float_t tempProb = 0.;
  entropyErr = TMath::Power(entropyErr,2.);
  Float_t tempErr = 0.;
  Float_t entropy = 0.;
  int *indices=new int[dim];
  for(int i=0;i<dim;i++) {
    indices[i]=1; //the first bin
  }
  while (1) {
    tempProb = aHist->GetBinContent(indices)/normWeight;
//    fprintf(stderr,"p(%i)=%.4f\t%.4f\n",indices[0],tempProb,aHist->GetBinContent(indices));
    tempErr = TMath::Power(aHist->GetBinError(indices)/normWeight,2) + TMath::Power(tempProb*normErr/normWeight,2);
    //tempErr(p) is sigma^2(p)
    if (tempProb) {
      entropy = entropy - tempProb*TMath::Log2(tempProb);
      entropyErr = entropyErr + TMath::Power(TMath::Abs(1+TMath::Log2(tempProb)),2)*tempErr;
      //entropyErr is sum_{bins} (1+log(p))^2*sigma^2(p)
      //changed this to reduce # of squares/sqrts
      }
    int j;
    for(j=0;j<dim;j++) {
      indices[j]++;
      if (indices[j] <= numBins[j]) break; //this is a valid index
      else indices[j]=1; //done with valid indices, loop around to 1 ( 0 is underflow)
    }
    if(j==dim) break;
  }
  entropyErr = sqrt(entropyErr);
  return entropy;

}

Float_t IntAndErr(THnSparseF*aHist, int dim, int *numBins,Float_t &integralErr) {
//  Double_t integral=aHist->ComputeIntegral(); // I don't know why this always returns 1
//  GetBinContent/Error works, so compute manually
  Float_t integral=0.;
  integralErr = 0;
  int *indices=new int[dim];
  for(int i=0;i<dim;i++) {
    indices[i]=1; //the first bin
  }
  while (1) {
    integral+=aHist->GetBinContent(indices);
    integralErr+=TMath::Power(aHist->GetBinError(indices),2.);
    int j;
    for(j=0;j<dim;j++) {
      indices[j]++;
      if (indices[j] <= numBins[j]) break; //this is a valid index
      else indices[j]=1; //done with valid indices, loop around
    }
    if(j==dim) break;
  }
 integralErr=sqrt(integralErr);
 return integral;   
}

/*void MutualInfoNVars(TString* theFormula = 0,TString histFileName = "fjet1MassSDb2", Int_t dim=1,Int_t* numBins = 0,Double_t* histMin = 0,Double_t* histMax = 0){
    if (dim==1) { // because vim syntax highlighting is stupid; don't ask
        if (theFormula==0) {
          theFormula = (TString*)malloc(1*sizeof(TString));
          theFormula[0]=TString("fjet1MassSDb2");
        }
        if (numBins==0) {
            numBins = (Int_t*)malloc(1*sizeof(Int_t));
            numBins[0]=80;
        }
        if (histMin==0) {
            histMin = (Double_t*)malloc(1*sizeof(Double_t));
            histMin[0]=-20;
        }
        if (histMax==0) {
            histMax = (Double_t*)malloc(1*sizeof(Double_t));
            histMax[0]=300;
        }
    }
*/
void MutualInfoNVars(TString cfgFileName="nVarConfig.txt") {
  //read config file
  Int_t dim=0;
  TString histFileName;
  ifstream cfgFile;
  cfgFile.open(cfgFileName.Data());
  cfgFile >> dim;
  if(!cfgFile.good()) { fprintf(stderr,"config file has bad format or does not exist %i\n",dim); exit(1); }
  cfgFile >> histFileName;
  TString*theFormula = new TString[dim];
  Int_t* numBins = new Int_t[dim];
  Double_t* histMin = new Double_t[dim];
  Double_t* histMax = new Double_t[dim];
  for(int i=0;i<dim;i++) { 
    cfgFile >> theFormula[i] >> numBins[i] >> histMin[i] >> histMax[i];
  }  
  cfgFile.close();

    //don't use so much memory you unmount hadoop
    Int_t expMem=1;
    for(int i=0;i<dim;i++)  expMem*=numBins[i];
    expMem*=6*sizeof(Float_t);
    fprintf(stderr,"WARNING: You are about to use up to %i B of memory\n",expMem);
    fprintf(stderr,"Pausing for %.2f seconds so you can re-evaluate and hit Ctrl-C",float(expMem)/(1024*1024*512));
    usleep(expMem*100/(1024*1024*512));

    //file io 
  TFile *signalFile = new TFile("signal_word.root");
  TTree *signalTree = (TTree*) signalFile->FindObjectAny("DMSTree");
  TFile *backgdFile = new TFile("background_word.root");
  TTree *backgdTree = (TTree*) backgdFile->FindObjectAny("DMSTree");
  TTreeFormula ** formulas = new TTreeFormula*[dim];
  for(int i=0;i<dim;i++) {    formulas[i] = new TTreeFormula(theFormula[i],theFormula[i],signalTree);  }
  TTreeFormula *fWeight = new TTreeFormula("weight","weight",signalTree);
  TTreeFormula *fSigCut = new TTreeFormula("abs(fjet1PartonId)","abs(fjet1PartonId)",signalTree);
    
  THnSparseF *signalHist = new THnSparseF("signalHist","signalHist",dim,numBins,histMin,histMax);  signalHist->Sumw2();
  THnSparseF *backgdHist = new THnSparseF("Backgd","Backgd",dim,numBins,histMin,histMax);  backgdHist->Sumw2();
  THnSparseF *sumHist = new THnSparseF("Sum","Sum",dim,numBins,histMin,histMax); sumHist->Sumw2();

  Double_t *evals = new Double_t[dim];
  Double_t weight;
  
  //fill signal
  int nEvents=signalTree->GetEntries();
  for(int i=0;i<nEvents;i++) {
    signalTree->GetEntry(i);
    if (fSigCut->EvalInstance()!=24) continue;
    for(int j=0;j<dim;j++) {
      evals[j] = formulas[j]->EvalInstance();
    }
    weight = fWeight->EvalInstance();
    signalHist->Fill(evals,weight);
    sumHist->Fill(evals,weight);
  }
  Float_t dummy;
  //fprintf(stderr,"%f %f\n",IntAndErr(signalHist,dim,numBins,dummy),test_int/IntAndErr(signalHist,dim,numBins,dummy));

  //background
  for(int i=0;i<dim;i++) {
    formulas[i]->Delete();
    formulas[i] = new TTreeFormula(theFormula[i],theFormula[i],backgdTree);
    //another broken thing in root: SetTree doesn't actually do anything
    //formulas[i]->SetTree(backgdTree);
  }
  fWeight->Delete();
  fWeight = new TTreeFormula("weight","weight",backgdTree);
  nEvents=backgdTree->GetEntries();
  for(int i=0;i<nEvents;i++) {
    backgdTree->GetEntry(i);
    for(int j=0;j<dim;j++) {
      evals[j] = formulas[j]->EvalInstance();
    }
    weight = fWeight->EvalInstance();
    backgdHist->Fill(evals,weight);
    sumHist->Fill(evals,weight);
  }

  //integrals and error prop
  Float_t signalErr = 0.;
  Float_t signalIntegral = IntAndErr(signalHist,1,numBins,signalErr);
  Float_t sumErr = 0.;
  Float_t sumIntegral = IntAndErr(sumHist,1,numBins,sumErr);
  Float_t signalFrac = signalIntegral/sumIntegral;
  Float_t signalFracErr = sqrt(TMath::Power(signalErr/sumIntegral,2) + TMath::Power(signalIntegral*sumErr/(TMath::Power(sumIntegral,2)),2));

  //entropies
  Float_t truthEntropy = -1*signalFrac*TMath::Log2(signalFrac) - (1-signalFrac)*TMath::Log2(1-signalFrac);
  Float_t truthErr = (TMath::Abs(1+TMath::Log2(signalFrac)) + TMath::Abs(1+TMath::Log2(1-signalFrac)))*signalFracErr;

  Float_t varErr = 0.;
  Float_t varEntropy = histEntropy(sumHist,dim,numBins,sumIntegral,sumErr,varErr);
  Float_t unionErr = 0.;
  Float_t unionEntropy = histEntropy(signalHist,dim,numBins,sumIntegral,sumErr,unionErr) +
                          histEntropy(backgdHist,dim,numBins,sumIntegral,sumErr,unionErr);
  //info
  Float_t MutualInfo = truthEntropy + varEntropy - unionEntropy;
  Float_t MutualErr = sqrt(TMath::Power(truthErr,2) + TMath::Power(varErr,2) + TMath::Power(unionErr,2));
  for(int i=0;i<dim;i++) {  cout << theFormula[i] << ((i==dim-1) ? endl : ":"); }
  cout << "H(T):   " << truthEntropy << " +- " << truthErr  << endl;
  cout << "H(A):   " << varEntropy   << " +- " << varErr    << endl;
  cout << "H(T,A): " << unionEntropy << " +- " << unionErr  << endl;
  cout << "I(T;A): " << MutualInfo   << " +- " << MutualErr << endl;

  //file io
  TFile *outFile = new TFile(TString("scratch/") + histFileName + TString(".root"),"RECREATE");
  signalHist->Write();
  backgdHist->Write();
  sumHist->Write();
  outFile->Close();

}
