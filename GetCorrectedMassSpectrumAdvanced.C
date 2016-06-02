//#include <iostream>
//using std::cout;
//using std::endl;

#include "TCanvas.h"
#include "TObject.h"
#include "TH1.h"
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

void GetCorrectedMassSpectrumAdvanced(TString fileName, TString histoOrCanvasName, Bool_t isMC){
  TFile *inputFile = new TFile(fileName,"READ");
  if ( !inputFile ) { cout<<"Fatal: the specified file can't be found."<<endl; return; }

  TObject *histoOrCanvas = dynamic_cast<TObject*>(inputFile->FindObjectAny(histoOrCanvasName.Data()));
  TString *className =new TString(((TClass*)histoOrCanvas->IsA())->GetName());

  if ( className->Contains("TH1") ){
    TH1 *histoToFit = dynamic_cast<TH1*>histoOrCanvas;
    if ( !histoToFit ) { cout<<"Fatal: the specified histogram can't be found."<<endl; return; }

    TF1* fitFunc;
    if ( isMC ) fitFunc = new TF1("fitFunc",FitFuncErfFixed,0.,12.,7);
    else fitFunc = new TF1("fitFunc",FitFuncErf,0.,12.,7);

    fitFunc->SetParameters(0.5, 1., 0.3, -1.e+03, 0.2, 0.1, 0.35);
    fitFunc->SetParNames("Norm", "p_{T} cut", "#sigma", "Offset", "Norm (LE)", "p_{T} cut (LE)", "#sigma (LE)", "Offset (LE)", "p_{T} change");

    histoToFit->Fit(fitFunc,"R");
    histoToFit->Draw("e");
  }

  TList *fittedHistos=new TList();

  if ( className->Contains("TCanvas") ){
    cout<<"Speicified object is a TCanvas"<<endl;
    TCanvas *canvasToFit = dynamic_cast<TCanvas*>histoOrCanvas;
    //canvasToFit->Draw();
    TList *listOfCanvasPrimitives = canvasToFit->GetListOfPrimitives();
    Int_t nCanvasPrimitives = listOfCanvasPrimitives->GetEntries();

    cout<<"canvas has "<<nCanvasPrimitives<<" primitives"<<endl;

    for (Int_t i = 0; i < nCanvasPrimitives; i++) {
      TObject *objectInCanvas = listOfCanvasPrimitives->At(i);

      if ( objectInCanvas->IsA() == TPad::Class() ) {
        cout<<"found pad: "<<i<<endl;
        TList *listOfPadPrimitives = dynamic_cast<TPad*>objectInCanvas)->GetListOfPrimitives();
        Int_t nPadPrimitives = listOfPadPrimitives->GetEntries();

        cout<<"pad has "<<nPadPrimitives<<" primitives"<<endl;

        TH1 *histoToFit;
        for (Int_t j = 0; j < nPadPrimitives; j++) {
          TObject *objectInPad = listOfPadPrimitives->At(j);

          if ( objectInPad->InheritsFrom(TH1::Class()) ){
            histoToFit = static_cast<TH1*>objectInPad;
            histoToFit->SetDirectory(0);
            cout<<"found histo: "<<histoToFit->GetName()<<endl;

            TF1* fitFunc = new TF1("fitFunc",FitFuncErfFixed,0.,12.,7);

            fitFunc->SetParameters(0.5, 1., 0.3, 1., 0.2, 0.1, 0.35);
            histoToFit->Fit(fitFunc,"RL");
            fittedHistos->Add(histoToFit);
            //histoToFit->DrawCopy("e");
            //canvasToFit->Update();
            //canvasToFit->Modified();



            histoToFit=0x0;
          }
        }
      }
    }
  }

  Int_t nBins = fittedHistos->GetEntries();

  Int_t canvasRows=0;
  Int_t canvasColumns=0;

  if ( nBins%2==0 ){
    canvasRows=2;
    canvasColumns=nBins/2;
  } else {
    if ( nBins%3==0 ){
      canvasRows=3;
      canvasColumns=nBins/3;
    } else {
      canvasRows=1;
      canvasColumns=nBins;
    }
  }

  TCanvas *canv=new TCanvas("canv","canv");
  canv->Divide(canvasColumns,canvasRows);

  for (Int_t i = 0; i < nBins; i++) {
    canv->cd(i+1);
    fittedHistos->At(i)->Draw();
  }
}
