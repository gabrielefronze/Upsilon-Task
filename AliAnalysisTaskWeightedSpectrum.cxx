#include <iostream>
using std::cout;
using std::endl;

#include "AliMCEvent.h"
#include "TH1.h"
#include "THnSparse.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TH1D.h"
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

AliAnalysisTaskWeightedSpectrum::AliAnalysisTaskWeightedSpectrum(Bool_t isMC, TString inputFileName, TString mode) :
  AliAnalysisTaskSE("analysis"),
  fOutput(0x0)

{
  fCuts=new AliMuonTrackCuts();
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;
  fIsMC=isMC;
  if ( mode.Contains("H") ) fMode="histo";
  if ( mode.Contains("F") ) fMode="function";

  cout<<"mode is "<<fMode<<endl;

  TFile *inputFile = new TFile(inputFileName,"READ");
  cout<<"file "<<inputFile<<endl;

  fRapidityAxis = new TAxis();
  inputFile->GetObject("rapidity_axis",fRapidityAxis);

  cout<<"got axis "<<fRapidityAxis->GetNbins()<<endl;

  fInputResponseFunctions =  new TObjArray();

  cout<<"created array "<<endl;

  TObject *buffer;

  if ( fMode=="histo" ){
    Double_t rapidityBinWidth = fRapidityAxis->GetBinWidth(1);
    Double_t rapidityLow = fRapidityAxis->GetXmin();

    cout<<"OK"<<endl;

    for (Int_t iInputResponseFunctions = 0; iInputResponseFunctions < fRapidityAxis->GetNbins(); iInputResponseFunctions++) {
      buffer = inputFile->FindObject(Form("Histo_ratio_MC_%f-%f",rapidityLow+rapidityBinWidth*iInputResponseFunctions,rapidityLow+rapidityBinWidth*iInputResponseFunctions+rapidityBinWidth));
      if (buffer) fInputResponseFunctions->Add(buffer);
      cout<<iInputResponseFunctions<<endl;
    }
  }

  if ( fMode=="function" ){
    for (Int_t iInputResponseFunctions = 0; iInputResponseFunctions < fRapidityAxis->GetNbins(); iInputResponseFunctions++) {
      buffer = inputFile->FindObject(Form("MCFittingFunction_%d",iInputResponseFunctions));
      if (buffer) fInputResponseFunctions->Add(buffer);
      cout<<iInputResponseFunctions<<endl;
    }
  }

  cout<<"Response parametrization obtained!"<<endl;

  buffer=0x0;
  //delete buffer;

  fTreeData = new Float_t[kWeight];

  inputFile->Close();

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskWeightedSpectrum::AliAnalysisTaskWeightedSpectrum(AliMuonTrackCuts *cuts, Bool_t isMC, TString inputFileName, TString mode) :
  AliAnalysisTaskSE("analysis"),
  fOutput(0x0)
{
  fCuts = cuts;
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;
  fIsMC=isMC;
  if ( mode.Contains("H") ) fMode="histo";
  if ( mode.Contains("F") ) fMode="function";

  cout<<"mode is "<<fMode<<endl;

  TFile *inputFile = new TFile(inputFileName,"READ");
  cout<<"file "<<inputFile<<endl;

  fRapidityAxis = new TAxis();
  inputFile->GetObject("rapidity_axis",fRapidityAxis);

  cout<<"got axis "<<fRapidityAxis->GetNbins()<<endl;

  fInputResponseFunctions =  new TObjArray();

  cout<<"created array "<<endl;

  TObject *buffer;

  if ( fMode=="histo" ){
    Double_t rapidityBinWidth = fRapidityAxis->GetBinWidth(1);
    Double_t rapidityLow = fRapidityAxis->GetXmin();

    cout<<"OK"<<endl;

    for (Int_t iInputResponseFunctions = 0; iInputResponseFunctions < fRapidityAxis->GetNbins(); iInputResponseFunctions++) {
      buffer = inputFile->FindObject(Form("Histo_ratio_MC_%f-%f",rapidityLow+rapidityBinWidth*iInputResponseFunctions,rapidityLow+rapidityBinWidth*iInputResponseFunctions+rapidityBinWidth));
      if (buffer) fInputResponseFunctions->Add(buffer);
      cout<<iInputResponseFunctions<<endl;
    }
  }

  if ( fMode=="function" ){
    for (Int_t iInputResponseFunctions = 0; iInputResponseFunctions < fRapidityAxis->GetNbins(); iInputResponseFunctions++) {
      buffer = inputFile->FindObject(Form("MCFittingFunction_%d",iInputResponseFunctions));
      if (buffer) fInputResponseFunctions->Add(buffer);
      cout<<iInputResponseFunctions<<endl;
    }
  }

  cout<<"Response parametrization obtained!"<<endl;

  buffer=0x0;
  //delete buffer;

  fTreeData = new Float_t[kWeight];

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
  if ( fInputResponseFunctions ) delete fInputResponseFunctions;
  if ( fTreeData ) delete fTreeData;
}


void AliAnalysisTaskWeightedSpectrum::NotifyRun()
{
  printf("Setting run number for cuts\n");
  fCuts->SetRun(fInputHandler);
}

void AliAnalysisTaskWeightedSpectrum::UserCreateOutputObjects()
{
  TTree *treeDimu=new TTree("MuonData","TTree containing analysis results");
  treeDimu->Branch("dimuon_momentum",&fTreeData[kMomentum]);
  treeDimu->Branch("dimuon_transverse_momentum",&fTreeData[kTransverse]);
  treeDimu->Branch("dimuon_rapidity",&fTreeData[kRapidity]);
  treeDimu->Branch("dimuon_mass",&fTreeData[kMass]);
  treeDimu->Branch("first_muon_transverse_momentum",&fTreeData[kPt1]);
  treeDimu->Branch("second_muon_transverse_momentum",&fTreeData[kPt2]);
  treeDimu->Branch("event_centrality",&fTreeData[kCentrality]);
  treeDimu->Branch("weight",&fTreeData[kWeight]);


  fOutput->AddAt(treeDimu,0);

  PostData(1,fOutput);
}

void AliAnalysisTaskWeightedSpectrum::UserExec(Option_t *)
{

  cout<<"ANALYSIS BEGINS"<<endl;
  if( !fIsMC ){
  	if( !(InputEvent()->GetFiredTriggerClasses()).Contains("CINT7-B-NOPF-MUFAST") ){
  		//cout<<"-> Rejected because of trigger class selection."<<endl;
  		return;
  	} else {
  		//cout<<"-> Accepted after applying trigger class selection"<<endl;
  	}
  }

  ////cout<<"Loading THnSparse"<<endl;
  TTree *treeDimu=(TTree*)fOutput->At(0);
  treeDimu->GetBranch("dimuon_momentum")->SetAddress(&fTreeData[kMomentum]);
  treeDimu->GetBranch("dimuon_transverse_momentum")->SetAddress(&fTreeData[kTransverse]);
  treeDimu->GetBranch("dimuon_rapidity")->SetAddress(&fTreeData[kRapidity]);
  treeDimu->GetBranch("dimuon_mass")->SetAddress(&fTreeData[kMass]);
  treeDimu->GetBranch("highest_muon_transverse_momentum")->SetAddress(&fTreeData[kPt1]);
  treeDimu->GetBranch("lowest_muon_transverse_momentum")->SetAddress(&fTreeData[kPt2]);
  treeDimu->GetBranch("event_centrality")->SetAddress(&fTreeData[kCentrality]);
  treeDimu->GetBranch("weight")->SetAddress(&fTreeData[kWeight]);

  // reading how much tracks are stored in the input event and the event centrality
  Int_t ntracks=AliAnalysisMuonUtility::GetNTracks(InputEvent());
  Double_t eventCentrality=0.;

  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  if ( multSelection ) eventCentrality=multSelection->GetMultiplicityPercentile("V0M");

  TObjArray *muPlus = new TObjArray();
  TObjArray *muMinus = new TObjArray();
  vector<Int_t> motherIndexes;

  // loop over the tracks to store only muon ones to obtain every possible dimuon
  AliAODTrack* muonBufferData=0x0;
  AliVParticle* muonBufferMC=0x0;
  AliVParticle* motherBufferMC=0x0;
  for(Int_t itrack=0;itrack<ntracks;itrack++){
    muonBufferData=(AliAODTrack*)AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());

// is the track  selected via standard muon cuts?
    if ( (!fCuts->IsSelected(muonBufferData))) continue;

    if ( AliAnalysisMuonUtility::GetMatchTrigger(muonBufferData)==0. ) continue;

    // il kFALSE è perchè in realtà la scelta dell'origine dei muoni servirà più tardi (quando valuteremo il sistematico)
    if( fIsMC && kFALSE ){ // if the analysed run is a MC run the macro excludes any particle not recognized as a muon from upsilon
      Int_t mcDaughterIndex=muonBufferData->GetLabel();
      Int_t mcMotherIndex=0;

      // check if there's any MC truth info related to the particle iTrack
      if ( mcDaughterIndex<0 ) continue; // this particle has no MC truth information
      else {
        muonBufferMC=MCEvent()->GetTrack(mcDaughterIndex);

        // check for particle identity
        if ( TMath::Abs(muonBufferMC->PdgCode())!=13 ) continue; // this particle is not a muon
        else {
          //cout<<"It's a muon! Cheers!!!"<<endl;
          mcMotherIndex=muonBufferMC->GetMother();
          if ( mcMotherIndex<0 ) continue;
          motherBufferMC=MCEvent()->GetTrack(mcMotherIndex);

          // check for particle mother identity
          if ( motherBufferMC->PdgCode()!=553 ) continue; // the mother of the studied particle is not a upsilon

          if ( muonBufferMC->Charge()>0. ) muPlus->AddAt(muonBufferData, mcMotherIndex);
          if ( muonBufferMC->Charge()<0. ) muMinus->AddAt(muonBufferData, mcMotherIndex);
          motherIndexes.push_back (mcMotherIndex);
          //else //cout<<"And its mother is an Upsilon! Double cheers!!!"<<endl;
        }
      }
    }
  }
  delete muonBufferData;
  delete muonBufferMC;
  delete motherBufferMC;

  AliAODTrack *firstMuon;
  AliAODTrack *secondMuon;
  for(std::vector<Int_t>::iterator mothersIterator = motherIndexes.begin(); mothersIterator != motherIndexes.end(); ++mothersIterator) {
    firstMuon=(AliAODTrack*)muPlus->At(*mothersIterator);
    secondMuon=(AliAODTrack*)muMinus->At(*mothersIterator);
    TLorentzVector dimuon;

    if (firstMuon->Charge()==secondMuon->Charge()) {
    	cout<<"Rejected (like-sign)"<<endl;
    	continue;
    }

    dimuon=AliAnalysisMuonUtility::GetTrackPair(firstMuon,secondMuon);
    fTreeData[kMomentum]=dimuon.P();
    fTreeData[kTransverse]=dimuon.Pt();
    fTreeData[kRapidity]=dimuon.Rapidity();
    fTreeData[kMass]=dimuon.M();
    fTreeData[kPt1]=firstMuon->Pt();
    fTreeData[kPt2]=secondMuon->Pt();
	  fTreeData[kCentrality]=(Double_t)eventCentrality;

    Int_t rapidityBin1 = fRapidityAxis->FindBin(firstMuon->Eta());
    Int_t rapidityBin2 = fRapidityAxis->FindBin(secondMuon->Eta());

    TObject *weightingFunction1 = fInputResponseFunctions->At(rapidityBin1);
    TObject *weightingFunction2 = fInputResponseFunctions->At(rapidityBin2);

    if ( fMode=="function" ) fTreeData[kWeight]=((TF1*)weightingFunction1)->Eval(fTreeData[kPt1])*((TF1*)weightingFunction2)->Eval(fTreeData[kPt2]);
    if ( fMode=="histo" ) fTreeData[kWeight]=((TH1D*)weightingFunction1)->GetBinContent(fTreeData[kPt1])*((TH1D*)weightingFunction2)->GetBinContent(fTreeData[kPt2]);

    //cout<<"Filling data"<<endl;
    if(fTreeData[kMass]>2.5)treeDimu->Fill();

   	fTreeData[kMomentum]=0.;
    fTreeData[kTransverse]=0.;
    fTreeData[kRapidity]=0.;
    fTreeData[kMass]=0.;
    fTreeData[kPt1]=0.;
    fTreeData[kPt2]=0.;
	  fTreeData[kCentrality]=0.;
    fTreeData[kWeight]=0.;
  }

  delete muonBufferData;
  delete muonBufferMC;
  delete motherBufferMC;

  PostData(1,fOutput);

  ////cout<<"########## ANALYSIS DONE! ##########"<<endl;

  return;
}

void AliAnalysisTaskWeightedSpectrum::Terminate(Option_t *) {
  fOutput=dynamic_cast<TList*>(GetOutputData(1));


  cout << "**********************" << endl;
  cout << "* Analysis completed *" << endl;
  cout << "**********************" << endl;

  return;
}
