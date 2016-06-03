#include <iostream>
using std::cout;
using std::endl;

#include "AliMCEvent.h"
#include "TH1.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TH3D.h"
#include "TList.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TObjArray.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliCentrality.h"
#include "AliAODTrack.h"
#include "TLorentzVector.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMultSelection.h"
#include "AliAnalysisTaskWeightedSpectrum.h"
#include "AliMuonTrackCuts.h"
#include <vector>

using namespace std;

ClassImp(AliAnalysisTaskWeightedSpectrum)

AliAnalysisTaskWeightedSpectrum::AliAnalysisTaskWeightedSpectrum() :
  AliAnalysisTaskSE("analysis"),
  fOutput(0x0),
  fCuts(NULL),
  fNEvents(0),
  fMode(""),
  fRapidityAxis(NULL),
  fInputResponseFunctionsMC(NULL),
  fInputResponseFunctionsData(NULL),
  fIsMC(0)
{
  cout<<"### default constructor"<<endl;
}


AliAnalysisTaskWeightedSpectrum::AliAnalysisTaskWeightedSpectrum(AliMuonTrackCuts *cuts, Bool_t isMC, TString inputFileName, TString mode) :
  AliAnalysisTaskSE("analysis"),
  fOutput(0x0)
{
  fOutput = new TList();
  fOutput->SetOwner();

  fCuts = cuts;
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;
  fIsMC=isMC;
  fMode=mode.Contains("F");

  cout<<"mode is "<<((fMode)?"function":"histo")<<endl;

  TFile *inputFile = new TFile(inputFileName,"READ");
  cout<<"file "<<inputFile<<endl;

  fRapidityAxis = new TAxis();
  inputFile->GetObject("rapidity_axis",fRapidityAxis);

  cout<<"got axis "<<fRapidityAxis->GetNbins()<<endl;

  fInputResponseFunctionsMC =  new TObjArray();
  fInputResponseFunctionsData =  new TObjArray();

  cout<<"created array "<<endl;

  TObject *buffer;

  if ( fMode==kFALSE ){
    Double_t rapidityBinWidth = fRapidityAxis->GetBinWidth(1);
    Double_t rapidityLow = fRapidityAxis->GetXmin();

    cout<<"OK"<<endl;

    for (Int_t iInputResponseFunctions = 0; iInputResponseFunctions < fRapidityAxis->GetNbins(); iInputResponseFunctions++) {

      printf("Histo_ratio_MC_%f-%f",rapidityLow+rapidityBinWidth*iInputResponseFunctions,rapidityLow+rapidityBinWidth*iInputResponseFunctions+rapidityBinWidth);

      inputFile->GetObject(Form("Histo_ratio_MC_%f-%f",rapidityLow+rapidityBinWidth*iInputResponseFunctions,rapidityLow+rapidityBinWidth*iInputResponseFunctions+rapidityBinWidth),buffer);
      if (!buffer) cout<<"problems"<<endl;
      static_cast<TH1D*>(buffer)->SetDirectory(0);
      fInputResponseFunctionsMC->Add(buffer);
      inputFile->GetObject(Form("Histo_ratio_data_%f-%f",rapidityLow+rapidityBinWidth*iInputResponseFunctions,rapidityLow+rapidityBinWidth*iInputResponseFunctions+rapidityBinWidth),buffer);
      if (!buffer) cout<<"problems"<<endl;
      static_cast<TH1D*>(buffer)->SetDirectory(0);
      fInputResponseFunctionsData->Add(buffer);
      cout<<iInputResponseFunctions<<endl;
    }
  }

  if ( fMode==kTRUE ){
    for (Int_t iInputResponseFunctions = 0; iInputResponseFunctions < fRapidityAxis->GetNbins(); iInputResponseFunctions++) {
      inputFile->GetObject(Form("MCFittingFunction_%d",iInputResponseFunctions),buffer);
      if (!buffer) cout<<"problems"<<endl;
      fInputResponseFunctionsMC->Add(buffer);
      inputFile->GetObject(Form("dataFittingFunction_%d",iInputResponseFunctions),buffer);
      if (!buffer) cout<<"problems"<<endl;
      fInputResponseFunctionsData->Add(buffer);
      cout<<iInputResponseFunctions<<endl;
    }
  }

  cout<<"Response parametrization obtained!"<<endl;

  buffer=0x0;
  //delete buffer;

  inputFile->Close();

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskWeightedSpectrum::~AliAnalysisTaskWeightedSpectrum()
{
  Info("~AliAnalysisTaskWeightedSpectrum","Calling Destructor");
  if ( fCuts ) fCuts=0x0;
  if ( fOutput ) delete fOutput;
  if ( fRapidityAxis ) delete fRapidityAxis;
  if ( fInputResponseFunctionsMC ) delete fInputResponseFunctionsMC;
  if ( fInputResponseFunctionsData ) delete fInputResponseFunctionsData;
}


void AliAnalysisTaskWeightedSpectrum::NotifyRun()
{
  printf("Setting run number for cuts\n");
  fCuts->SetRun(fInputHandler);
}

void AliAnalysisTaskWeightedSpectrum::UserCreateOutputObjects()
{

  cout<<"Creating Output objects"<<endl;
  TH3D *histoMC = new TH3D("histo_MC","histo_MC",1000,0.,12.,20,2.5,4.,1000,0.,50.);
  TH3D *histoData = new TH3D("histo_data","histo_data",1000,0.,12.,20,2.5,4.,1000,0.,50.);


  fOutput->AddAt(histoMC,0);
  fOutput->AddAt(histoData,1);

  cout<<"Output correctly instanced"<<endl;

  PostData(1,fOutput);
}

void AliAnalysisTaskWeightedSpectrum::UserExec(Option_t *)
{

  //cout<<"ANALYSIS BEGINS"<<endl;
  TH3D *histoMC=static_cast<TH3D*>(fOutput->At(0));
  TH3D *histoData=static_cast<TH3D*>(fOutput->At(1));
  //cout<<"Loading THnSparse"<<endl;


  // reading how much tracks are stored in the input event and the event centrality
  Int_t ntracks=AliAnalysisMuonUtility::GetNTracks(InputEvent());
  Double_t eventCentrality=0.;

  //cout<<"Centrality computation..."<<endl;
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  if ( multSelection ) eventCentrality=multSelection->GetMultiplicityPercentile("V0M");

  TObjArray *muPlus = new TObjArray();
  TObjArray *muMinus = new TObjArray();
  vector<Int_t> motherIndexes;

  //cout<<"Created TObjArrays"<<endl;

  // loop over the tracks to store only muon ones to obtain every possible dimuon
  AliAODTrack* muonBufferData=0x0;
  AliVParticle* muonBufferMC=0x0;
  AliVParticle* motherBufferMC=0x0;
  for(Int_t itrack=0;itrack<ntracks;itrack++){

    //cout<<"### "<<itrack<<endl;

    muonBufferData=(AliAODTrack*)AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());

    // is the track  selected via standard muon cuts?
    if ( (!fCuts->IsSelected(muonBufferData))) continue;

    // il kFALSE è perchè in realtà la scelta dell'origine dei muoni servirà più tardi (quando valuteremo il sistematico)
    if( fIsMC ){ // if the analysed run is a MC run the macro excludes any particle not recognized as a muon from upsilon
      Int_t mcDaughterIndex=muonBufferData->GetLabel();
      Int_t mcMotherIndex=0;

      // check if there's any MC truth info related to the particle iTrack
      if ( mcDaughterIndex<0 ) continue; // this particle has no MC truth information
      else {
        muonBufferMC=MCEvent()->GetTrack(mcDaughterIndex);

        // check for particle identity
        if ( TMath::Abs(muonBufferMC->PdgCode())!=13 ) continue; // this particle is not a muon
        //cout<<"It's a muon! Cheers!!!"<<endl;

        mcMotherIndex=muonBufferMC->GetMother();
        if ( mcMotherIndex<0 ) continue;
        motherBufferMC=MCEvent()->GetTrack(mcMotherIndex);

        // check for particle mother identity
        if ( motherBufferMC->PdgCode()!=553 ) continue; // the mother of the studied particle is not a upsilon
        if ( muonBufferMC->Charge()>0. ) muPlus->AddAt(muonBufferData, mcMotherIndex);
        if ( muonBufferMC->Charge()<0. ) muMinus->AddAt(muonBufferData, mcMotherIndex);
        motherIndexes.push_back (mcMotherIndex);
        //cout<<"found upsilon #"<<mcMotherIndex<<endl;
        //else //cout<<"And its mother is an Upsilon! Double cheers!!!"<<endl;
      }
    }
  }

  AliAODTrack *firstMuon;
  AliAODTrack *secondMuon;
  for(std::vector<Int_t>::iterator mothersIterator = motherIndexes.begin(); mothersIterator != motherIndexes.end(); ++mothersIterator) {
    //cout<<"mother #"<<*mothersIterator<<endl;
    firstMuon=(AliAODTrack*)muPlus->At(*mothersIterator);
    secondMuon=(AliAODTrack*)muMinus->At(*mothersIterator);

    if ( !firstMuon || !secondMuon ) continue;

    TLorentzVector dimuon = AliAnalysisMuonUtility::GetTrackPair(firstMuon,secondMuon);
    //cout<<"dimuon created"<<endl;
    Int_t rapidityBin1 = fRapidityAxis->FindBin(firstMuon->Eta());
    Int_t rapidityBin2 = fRapidityAxis->FindBin(secondMuon->Eta());
    Double_t pt1 = firstMuon->Pt();
    Double_t pt2 = secondMuon->Pt();

    TObject *weightingFunction1MC = fInputResponseFunctionsMC->At(rapidityBin1);
    TObject *weightingFunction2MC = fInputResponseFunctionsMC->At(rapidityBin2);

    TObject *weightingFunction1Data = fInputResponseFunctionsData->At(rapidityBin1);
    TObject *weightingFunction2Data = fInputResponseFunctionsData->At(rapidityBin2);

    Double_t weightMC;
    Double_t weightData;

    //cout<<"Computing weights"<<endl;
    //cout<<fInputResponseFunctionsMC->GetEntries()<<endl;

    if ( fMode==kFALSE ){

      //cout<<"histo mode"<<endl;
      TH1D *weightingHisto1MC = static_cast<TH1D*>(weightingFunction1MC);
      TH1D *weightingHisto1Data = static_cast<TH1D*>(weightingFunction1Data);
      TH1D *weightingHisto2MC = static_cast<TH1D*>(weightingFunction2MC);
      TH1D *weightingHisto2Data = static_cast<TH1D*>(weightingFunction2Data);

      weightMC = weightingHisto1MC->GetBinContent(weightingHisto1MC->FindBin(pt1));
      weightData = weightingHisto1MC->GetBinContent(weightingHisto1Data->FindBin(pt1));
      weightMC *= weightingHisto2MC->GetBinContent(weightingHisto2MC->FindBin(pt2));
      weightData *= weightingHisto2MC->GetBinContent(weightingHisto2Data->FindBin(pt2));

    } else if ( fMode==kTRUE ){

      //cout<<"function mode"<<endl;

      TF1 *weightingFunctionMC1 = static_cast<TF1*>(weightingFunction1MC);
      TF1 *weightingFunctionData1 = static_cast<TF1*>(weightingFunction1Data);
      TF1 *weightingFunctionMC2 = static_cast<TF1*>(weightingFunction2MC);
      TF1 *weightingFunctionData2 = static_cast<TF1*>(weightingFunction2Data);

      weightMC = weightingFunctionMC1->Eval(pt1);
      weightData = weightingFunctionData1->Eval(pt1);

      weightMC *= weightingFunctionMC2->Eval(pt2);
      weightData *= weightingFunctionData2->Eval(pt2);
    }

    histoMC->Fill(dimuon.Pt(), dimuon.Rapidity(), dimuon.M(),weightMC);
    histoData->Fill(dimuon.Pt(), dimuon.Rapidity(), dimuon.M(),weightData);
    //cout<<"filled"<<endl;

  }

  delete muPlus;
  delete muMinus;

  PostData(1,fOutput);

  ////cout<<"########## ANALYSIS DONE! ##########"<<endl;

  return;
}

void AliAnalysisTaskWeightedSpectrum::Terminate(Option_t *) {
  fOutput=dynamic_cast<TList*>(GetOutputData(1));

  TH3D *histoMC=static_cast<TH3D*>(fOutput->At(0));
  TH3D *histoData=static_cast<TH3D*>(fOutput->At(1));

  TCanvas *canv=new TCanvas("canv","canv");
  canv->Divide(1,2);
  canv->cd(1);
  histoMC->Draw("LEGO2Z");
  canv->cd(2);
  histoData->Draw("LEGO2Z");


  cout << "**********************" << endl;
  cout << "* Analysis completed *" << endl;
  cout << "**********************" << endl;

  return;
}
