#include <iostream>
using std::cout;
using std::endl;

//#include "TH1.h"
#include "THnSparse.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TH1D.h"
#include "TF1.h"
#include "TList.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TObjArray.h"

//______________________________________________
Double_t FitFuncErfFixed ( Double_t* xVal, Double_t* par )
{
  Double_t xx = xVal[0];
  Double_t currX = TMath::Max(xx,par[6]);
  Double_t sqrtTwo = TMath::Sqrt(2.);
  Double_t yVal = 1.+par[0]*(TMath::Erf((currX-par[1])/par[2]/sqrtTwo)-1.);
  if ( xx < par[6] ) yVal += par[3]*(TMath::Erf((-xx-par[4])/par[5]/sqrtTwo) - TMath::Erf((-par[6]-par[4])/par[5]/sqrtTwo));

  return yVal;
}

//______________________________________________
Double_t FitFuncErf ( Double_t* xVal, Double_t* par )
{
  Double_t xx = xVal[0];
  Double_t yVal = 0.;
  if ( xx > par[8] )
    yVal = par[3]+par[0]*(1.+TMath::Erf((xx-par[1])/par[2]/TMath::Sqrt(2.)));
  else
    yVal = par[7]+par[4]*(1.+TMath::Erf((-xx-par[5])/par[6]/TMath::Sqrt(2.)));

  return yVal;
}

void GetCorrectedMassSpectrum(TString fileName, TString histoName, Bool_t isMC){
  TFile *inputFile = new TFile(Form("UpsilonTaskOutputRatiosSparse/ResponseFunctions%s.root",((isMC)?"_MC":"_data")),"READ");
  if ( !inputFile ) { cout<<"Fatal: the specified file can't be found."<<endl; return; }

  TH1D *histoToFit = new TH1D();
  histoToFit = dynamic_cast<TH1D*>(inputFile->FindObjectAny(histoName.Data()));
  if ( !histoToFit ) { cout<<"Fatal: the specified histogram can't be found."<<endl; return; }

  TF1* fitFunc;
  if ( isMC ) fitFunc = new TF1("fitFunc",FitFuncErfFixed,0.,10.,7);
  else fitFunc = new TF1("fitFunc",FitFuncErf,0.,10.,7);

  fitFunc->SetParameters(0.5, 1., 0.3, 1., 0.2, 0.1, 0.35);

  histoToFit->Fit(fitFunc,"R");
  histoToFit->Draw("e");

}
