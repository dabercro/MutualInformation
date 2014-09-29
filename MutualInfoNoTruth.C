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
#include "THnSparse.h"
#include "TTreeFormula.h"

//identical to MutualInfoNVars, but computes I(a,b) instead of I(T,A)
//where a,b are single variables instead of sets of variables
//to be used for benchmarking mutual info with truth

Float_t histEntropy(THnSparseF *aHist,int dim, int* numBins,Float_t normWeight,Float_t normErr,Float_t& entropyErr){

  Float_t tempProb = 0.;
  entropyErr = TMath::Power(entropyErr,2.);
  Float_t tempErr = 0.;
  Float_t entropy = 0.;
  Double_t minusInfinity=TMath::Log2(0);
  int *indices=new int[dim];
  for(int i=0;i<dim;i++) {
    indices[i]=1; //the first bin
  }
  while (1) {
    tempProb = aHist->GetBinContent(indices)/normWeight;
    tempErr = TMath::Power(aHist->GetBinError(indices)/normWeight,2) + TMath::Power(tempProb*normErr/normWeight,2);
    //tempErr(p) is sigma^2(p)
    Double_t logProb = TMath::Log2(tempProb);
    if (logProb>minusInfinity) {
      entropy = entropy - tempProb*logProb;
      entropyErr = entropyErr + TMath::Power(TMath::Abs(1+logProb),2)*tempErr;
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

void MutualInfoNoTruth(TString cfgFileName="noTruthConfig.txt") {
  //read config file
  Int_t dim;
  TString histFileName;
  ifstream cfgFile;
  cfgFile.open(cfgFileName.Data());
  cfgFile >> histFileName;
  cfgFile >> dim;
  if (dim!=2) { fprintf(stderr,"this script should only be used to compare two variables at a time\n"); exit(1); }
  if(!cfgFile.good()) { fprintf(stderr,"config file has bad format or does not exist %i\n",dim); exit(1); }
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
    fprintf(stderr,"WARNING: You are about to use up to %.2f MB of memory\n",float(expMem)/(1024*1024));
    fprintf(stderr,"Pausing for %.2f seconds so you can re-evaluate and hit Ctrl-C",float(expMem)/(1024*1024*512));
    usleep(expMem*100/(1024*1024*512));
    cout << endl << endl;

    //file io 
  TFile *signalFile = new TFile("signal_word.root");
  TTree *signalTree = (TTree*) signalFile->FindObjectAny("DMSTree");
  TFile *backgdFile = new TFile("background_word.root");
  TTree *backgdTree = (TTree*) backgdFile->FindObjectAny("DMSTree");
  TTreeFormula ** formulas = new TTreeFormula*[dim];
  for(int i=0;i<dim;i++) {    formulas[i] = new TTreeFormula(theFormula[i],theFormula[i],signalTree);  }
  TTreeFormula *fWeight = new TTreeFormula("weight","weight",signalTree);
  TTreeFormula *fSigCut = new TTreeFormula("abs(fjet1PartonId)","abs(fjet1PartonId)",signalTree);
    
  THnSparseF *AHist = new THnSparseF("AHist","AHist",1,&(numBins[0]),&(histMin[0]),&(histMax[0]));  AHist->Sumw2();
  THnSparseF *BHist = new THnSparseF("BHist","BHist",1,&(numBins[1]),&(histMin[1]),&(histMax[1]));  BHist->Sumw2();
  THnSparseF *ABHist = new THnSparseF("ABHist","ABHist",2,numBins,histMin,histMax); ABHist->Sumw2();

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
    AHist->Fill(&(evals[0]),weight);
    BHist->Fill(&(evals[1]),weight);
    ABHist->Fill(evals,weight);
  }
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
    AHist->Fill(&(evals[0]),weight);
    BHist->Fill(&(evals[1]),weight);
    ABHist->Fill(evals,weight);
  }

  //integrals and error prop
  Float_t AIntErr = 0.;
  Float_t AIntegral = IntAndErr(AHist,1,&(numBins[0]),AIntErr);
  Float_t BIntErr = 0.;
  Float_t BIntegral = IntAndErr(BHist,1,&(numBins[1]),BIntErr);
  Float_t ABIntErr = 0.;
  Float_t ABIntegral = IntAndErr(ABHist,dim,numBins,ABIntErr);
  //entropies

  Float_t AErr = 0.;
  Float_t AEntropy = histEntropy(AHist,1,&(numBins[0]),AIntegral,AIntErr,AErr);
  Float_t BErr = 0.;
  Float_t BEntropy = histEntropy(BHist,1,&(numBins[1]),BIntegral,BIntErr,BErr);
  Float_t ABErr = 0.;
  Float_t ABEntropy = histEntropy(ABHist,dim,&(numBins[1]),ABIntegral,ABIntErr,ABErr);

  //info
  Float_t MutualInfo = AEntropy+BEntropy-ABEntropy;
  Float_t MutualErr = sqrt(TMath::Power(AErr,2) + TMath::Power(BErr,2) + TMath::Power(ABErr,2));
  for(int i=0;i<dim;i++) {  
    cout << theFormula[i] ;
    if (i==dim-1) cout<< endl; 
    else cout<< ":";
  }
  cout << "H(A):   " << AEntropy << " +- " << AErr  << endl;
  cout << "H(B):   " << BEntropy   << " +- " << BErr    << endl;
  cout << "H(A,B): " << ABEntropy << " +- " << ABErr  << endl;
  cout << "I(A;B): " << MutualInfo   << " +- " << MutualErr << endl;

  //file io
  TFile *outFile = new TFile(TString("scratch/") + histFileName + TString(".root"),"RECREATE");
  AHist->Write();
  BHist->Write();
  ABHist->Write();
  outFile->Close();

}
