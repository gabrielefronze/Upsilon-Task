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

Bool_t CheckTObjectPointer(TObject *obj){
  if ( !obj ) {
    printf("FATAL: Pointer named %s of class %s is not valid.",obj->GetName(),obj->IsA()->GetName());
    return kTRUE;
  }
  return kFALSE;
}

void CreateHistrogramFromSparse(TString analysisResultPath, TString sparseName){

  TFile *analysisResult = new TFile(Form("%s/AnalysisResults.root",analysisResultPath.Data()),"READ");
  TFile *outputFile = new TFile(Form("%s/AptOverLptHistograms.root",analysisResultPath.Data()),"RECREATE");
  if ( CheckTObjectPointer(analysisResult) ) return;
  if ( CheckTObjectPointer(outputFile) ) return;

  THnSparse *dataSparse = static_cast<THnSparse*>(analysisResult->FindObject(sparseName));
  if ( CheckTObjectPointer(dataSparse) ) return;

  Int_t numberOfRapidityBins=10;
  Double_t rapidityLow=2.5;
  Double_t rapidityHi=4.0;
  Double_t rapidityStep=(rapidityHi-rapidityLow)/numberOfRapidityBins;

}
