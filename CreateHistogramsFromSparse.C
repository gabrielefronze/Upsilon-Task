#include <iostream>
#include "THnSparse.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TList.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TObjArray.h"

using namespace std;

enum{
  kTriggerFlag,
  kPt,
  kRapidity,
  kCentrality,
  kLocalBoard,
  kSparseDimension
};

Bool_t CheckTObjectPointer(TObject *obj){
  if ( !obj ) {
    printf("FATAL: Pointer named %s of class %s is not valid.",obj->GetName(),obj->IsA()->GetName());
    return kTRUE;
  }
  return kFALSE;
}

Double_t FitFuncErfFixed ( Double_t* xVal, Double_t* par )
{
  Double_t xx = xVal[0];
  Double_t currX = TMath::Max(xx,par[6]);
  Double_t sqrtTwo = TMath::Sqrt(2.);
  Double_t yVal = 1.+par[0]*(TMath::Erf((currX-par[1])/par[2]/sqrtTwo)-1.);
  if ( xx < par[6] ) yVal += par[3]*(TMath::Erf((-xx-par[4])/par[5]/sqrtTwo) - TMath::Erf((-par[6]-par[4])/par[5]/sqrtTwo));

  return yVal;
}

void CreateHistogramsFromSparse(TString dataAnalysisResultPath, TString dataSparseName, TString MCAnalysisResultPath, TString MCSparseName){

  printf("Opening files... Output can be found in %s/AptOverLptHistograms.root\n",dataAnalysisResultPath.Data());
  TFile *dataAnalysisResult = new TFile(Form("%s/AnalysisResults.root",dataAnalysisResultPath.Data()),"READ");
  TFile *MCAnalysisResult = new TFile(Form("%s/AnalysisResults.root",MCAnalysisResultPath.Data()),"READ");
  TFile *outputFile = new TFile(Form("%s/AptOverLptHistograms.root",dataAnalysisResultPath.Data()),"RECREATE");
  if ( CheckTObjectPointer(dataAnalysisResult) ) return;
  if ( CheckTObjectPointer(outputFile) ) return;

  printf("Retrieving THnSparses..\n");
  THnSparse *dataSparse = static_cast<THnSparse*>(dataAnalysisResult->FindObject(dataSparseName));
  if ( CheckTObjectPointer(dataSparse) ) return;
  THnSparse *MCSparse = static_cast<THnSparse*>(MCAnalysisResult->FindObject(MCSparseName));
  if ( CheckTObjectPointer(MCSparse) ) return;

  const Int_t numberOfRapidityBins=10;
  Int_t rebinFactor=10;

  Double_t rapidityLow=2.5;
  Double_t rapidityHi=4.0;
  Double_t binWidth=(rapidityHi-rapidityLow)/numberOfRapidityBins;

  printf("Analysis will be performed with %d bins in rapidity and a rebin of %d\n",numberOfRapidityBins,rebinFactor);

  TH1D *histosAptBuffer;

  TH1D *dataHistosRatio[numberOfRapidityBins];
  TF1 *dataFittingFunctions[numberOfRapidityBins];

  printf("Creating histograms from real data\n");
  for (Int_t i = 0; i < numberOfRapidityBins; i++) {
    Double_t binLow=rapidityLow+i*binWidth;
    dataSparse->GetAxis(kRapidity)->SetRangeUser(binLow, binLow+binWidth);

    dataSparse->GetAxis(kTriggerFlag)->SetRangeUser(1., 1.);
    histosAptBuffer=dataSparse->Projection(kPt,"E"); // this is the Apt distribution
    histosAptBuffer->Rebin(rebinFactor);

    dataSparse->GetAxis(kTriggerFlag)->SetRangeUser(2., 2.);
    dataHistosRatio[i]=dataSparse->Projection(kPt,"E"); // this is the Lpt distribution but will be divided by Apt
    dataHistosRatio[i]->Rebin(rebinFactor);

    (dataHistosRatio[i])->SetName(Form("Histo_ratio_data_%f-%f",binLow,binLow+binWidth));
    (dataHistosRatio[i])->SetTitle(Form("Histo ratio data %f-%f",binLow,binLow+binWidth));
    (dataHistosRatio[i])->Divide(histosAptBuffer); // now this histogram contains the ratio of the two
    (dataHistosRatio[i])->SetLineColor(kRed);

    delete histosAptBuffer;

    dataFittingFunctions[i] = new TF1(Form("dataFittingFunction_%d",i),FitFuncErfFixed,0.,10.,7);
    dataFittingFunctions[i]->SetParameters(0.5, 1., 0.3, 1., 0.2, 0.1, 0.35);
    (dataHistosRatio[i])->Fit(dataFittingFunctions[i],"RL");
  }


  TH1D *MCHistosRatio[numberOfRapidityBins];

  printf("Creating histograms from MC data\n");
  for (Int_t i = 0; i < numberOfRapidityBins; i++) {
    Double_t binLow=rapidityLow+i*binWidth;
    MCSparse->GetAxis(kRapidity)->SetRangeUser(binLow, binLow+binWidth);

    MCSparse->GetAxis(kTriggerFlag)->SetRangeUser(1., 1.);
    histosAptBuffer=MCSparse->Projection(kPt,"E"); // this is the Apt distribution
    histosAptBuffer->Rebin(rebinFactor);

    MCSparse->GetAxis(kTriggerFlag)->SetRangeUser(2., 2.);
    MCHistosRatio[i]=MCSparse->Projection(kPt,"E"); // this is the Lpt distribution but will be divided by Apt
    MCHistosRatio[i]->Rebin(rebinFactor);

    (MCHistosRatio[i])->SetName(Form("Histo_ratio_MC_%f-%f",binLow,binLow+binWidth));
    (MCHistosRatio[i])->SetTitle(Form("Histo ratio MC %f-%f",binLow,binLow+binWidth));
    (MCHistosRatio[i])->Divide(histosAptBuffer); // now this histogram contains the ratio of the two
    (MCHistosRatio[i])->SetLineColor(kBlack);

    delete histosAptBuffer;

    TF1* MCFittingFunction = new TF1(Form("dataFittingFunction_%d",i),FitFuncErfFixed,0.,10.,7);
    MCFittingFunction->SetParameters(0.5, 1., 0.3, 1., 0.2, 0.1, 0.35);
    MCFittingFunction->FixParameter(7, dataFittingFunctions[i]->GetParameter(7));

    // here the parameter obtained from data has to be fixed!
    (MCHistosRatio[i])->Fit(dataFittingFunctions[i],"RL");
  }

  Int_t canvasRows=0;
  Int_t canvasColumns=0;

  if ( numberOfRapidityBins%2==0 ){
    canvasRows=2;
    canvasColumns=numberOfRapidityBins/2;
  } else {
    if ( numberOfRapidityBins%3==0 ){
      canvasRows=3;
      canvasColumns=numberOfRapidityBins/3;
    } else {
      canvasRows=1;
      canvasColumns=numberOfRapidityBins;
    }
  }

  TCanvas *canvasRatios = new TCanvas("Lpt_over_Apt_ratios","Lpt over Apt ratios");
  canvasRatios->Divide(canvasColumns,canvasRows);

  for (Int_t iHistos = 0; iHistos < numberOfRapidityBins; iHistos++) {
    canvasRatios->cd(iHistos+1);
    (dataHistosRatio[iHistos])->Draw("E");
    (MCHistosRatio[iHistos])->Draw("SAME E");

    outputFile->cd();
    (dataHistosRatio[iHistos])->Write();
    (MCHistosRatio[iHistos])->Write();
  }

}
