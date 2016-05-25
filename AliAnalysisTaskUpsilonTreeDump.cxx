#include <iostream>
using std::cout;
using std::endl;

#include "TH1.h"
#include "TMath.h"
#include "TTree.h"
#include "TList.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TObjArray.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliCentrality.h"
#include "AliAODTrack.h"
#include "TLorentzVector.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMultSelection.h"
#include "AliAnalysisTaskUpsilonTreeDump.h"
#include "AliMuonTrackCuts.h"


ClassImp(AliAnalysisTaskUpsilonTreeDump
)

AliAnalysisTaskUpsilonTreeDump::AliAnalysisTaskUpsilonTreeDump() :
  AliAnalysisTaskSE(""), 
  fAODEvent(0x0),
  fESDEvent(0x0),
  fOutput(0x0)
{
  fCuts=new AliMuonTrackCuts();
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;
  fTreeData=(Float_t*)malloc(kSparseDimension*sizeof(Float_t));
  fPassedCutsMuon=(Bool_t*)malloc(kDimension*sizeof(Bool_t));

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskUpsilonTreeDump::AliAnalysisTaskUpsilonTreeDump(const char *name, AliMuonTrackCuts *cuts) :
  AliAnalysisTaskSE(name),
  fAODEvent(0x0),
  fESDEvent(0x0),
  fOutput(0x0)
{
  fCuts=cuts;
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;
  fTreeData=(Float_t*)malloc(kSparseDimension*sizeof(Float_t));
  fPassedCutsMuon=(Bool_t*)malloc(kDimension*sizeof(Bool_t));

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskUpsilonTreeDump::~AliAnalysisTaskUpsilonTreeDump()
{
  Info("~AliAnalysisTaskUpsilonTreeDump","Calling Destructor");
  if(fCuts)fCuts=0x0;
  if(fOutput)delete fOutput;
  if(fTreeData)delete fTreeData;
  if(fPassedCutsMuon)delete fPassedCutsMuon;
}


void AliAnalysisTaskUpsilonTreeDump::NotifyRun()
{
  printf("Setting run number for cuts\n");
  fCuts->SetRun(fInputHandler);
}

void AliAnalysisTaskUpsilonTreeDump::UserCreateOutputObjects()
{
  fOutput = new TList();
  fOutput->SetOwner();

  TTree *treeMu=new TTree("MuonData","TTree containing analysis results");
  treeMu->Branch("muon_momentum",&fTreeData[kMomentum]);
  treeMu->Branch("muon_transverse_momentum",&fTreeData[kTransverse]);
  treeMu->Branch("muon_rapidity",&fTreeData[kRapidity]);
  treeMu->Branch("muon_charge",&fTreeData[kCharge]);
  treeMu->Branch("event_centrality",&fTreeData[kCentrality]);

  treeMu->Branch("Physics_Selection",&fPassedCutsMuon[kPhysicsSelection]);
  treeMu->Branch("Trigger_Selection",&fPassedCutsMuon[kTriggerSelection]);
  treeMu->Branch("Has_Tracker_Data",&fPassedCutsMuon[kHasTrackerData]);
  treeMu->Branch("Low_Pt_Cut",&fPassedCutsMuon[kLowPtCut]);
  treeMu->Branch("Eta_Cut",&fPassedCutsMuon[kEtaCut]);
  treeMu->Branch("RAbs_Cut",&fPassedCutsMuon[kRAbsCut]);
  treeMu->Branch("Standard_Muon_Cuts",&fPassedCutsMuon[kStandardMuonCuts]);

  treeMu->Branch("Event_Number",&fNEvents);

  fOutput->AddAt(treeMu,0);

  PostData(1,fOutput);
} 


void AliAnalysisTaskUpsilonTreeDump::UserExec(Option_t *)
{  
  // sparse that contains data about the muon data
  TTree* treeMu=((TTree*)fOutput->At(0));
  treeMu->GetBranch("muon_momentum")->SetAddress(&fTreeData[kMomentum]);
  treeMu->GetBranch("muon_transverse_momentum")->SetAddress(&fTreeData[kTransverse]);
  treeMu->GetBranch("muon_rapidity")->SetAddress(&fTreeData[kRapidity]);
  treeMu->GetBranch("muon_charge")->SetAddress(&fTreeData[kCharge]);
  treeMu->GetBranch("event_centrality")->SetAddress(&fTreeData[kCentrality]);

  treeMu->GetBranch("Physics_Selection")->SetAddress(&fPassedCutsMuon[kPhysicsSelection]);
  treeMu->GetBranch("Trigger_Selection")->SetAddress(&fPassedCutsMuon[kTriggerSelection]);
  treeMu->GetBranch("Has_Tracker_Data")->SetAddress(&fPassedCutsMuon[kHasTrackerData]);
  treeMu->GetBranch("Low_Pt_Cut")->SetAddress(&fPassedCutsMuon[kLowPtCut]);
  treeMu->GetBranch("Eta_Cut")->SetAddress(&fPassedCutsMuon[kEtaCut]);
  treeMu->GetBranch("RAbs_Cut")->SetAddress(&fPassedCutsMuon[kRAbsCut]);
  treeMu->GetBranch("Standard_Muon_Cuts")->SetAddress(&fPassedCutsMuon[kStandardMuonCuts]);

  treeMu->GetBranch("Event_Number")->SetAddress(&fNEvents);

	cout<<"Run:"<<InputEvent()->GetRunNumber()<<" Event:"<<fNEvents++;
  Int_t treeEntries=treeMu->GetEntries();

  for(Int_t iCuts=0; iCuts<kDimension; iCuts++){
    fPassedCutsMuon[iCuts]=kFALSE;
  }

  fPassedCutsMuon[kPhysicsSelection]=kTRUE;
  fPassedCutsMuon[kTriggerSelection]=((InputEvent()->GetFiredTriggerClasses()).Contains("CMUL7-B-NOPF-MUFAST"));
  
  // reading how much tracks are stored in the input event and the event centrality
  Int_t ntracks=AliAnalysisMuonUtility::GetNTracks(InputEvent());
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  fTreeData[kCentrality]=multSelection->GetMultiplicityPercentile("V0M");

  // loop over the tracks to store only muon ones to obtain every possible muon
  AliAODTrack* muonBuffer=0x0;
  for(Int_t itrack=0;itrack<ntracks;itrack++){
    muonBuffer=(AliAODTrack*)AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());

    // check if the track is seen in the muon tracker
    fPassedCutsMuon[kHasTrackerData]=AliAnalysisMuonUtility::IsMuonTrack(muonBuffer);

    // keep only low-pt matching tracks
    fPassedCutsMuon[kLowPtCut]=AliAnalysisMuonUtility::GetMatchTrigger(muonBuffer)>2;

    // Rapidity acceptance cut.
    fPassedCutsMuon[kEtaCut]=(muonBuffer->Eta()>-4.0 && muonBuffer->Eta()<-2.5);

    // Rabs cut.
    fPassedCutsMuon[kRAbsCut]=(muonBuffer->GetRAtAbsorberEnd()>17.6 && muonBuffer->GetRAtAbsorberEnd()<89.5);

    // is the track  selected via standard muon cuts?
    fPassedCutsMuon[kStandardMuonCuts]=fCuts->IsSelected(muonBuffer);

    fTreeData[kMomentum]=muonBuffer->P();
    fTreeData[kTransverse]=muonBuffer->Pt();
    fTreeData[kRapidity]=muonBuffer->Eta();
    fTreeData[kCharge]=muonBuffer->Charge();

    treeMu->Fill();

    for(Int_t iData=0; iData<kSparseDimension; iData++){
      fTreeData[iData]=0.;
    }

    fPassedCutsMuon[kHasTrackerData]=kFALSE;
    fPassedCutsMuon[kLowPtCut]=kFALSE;
    fPassedCutsMuon[kEtaCut]=kFALSE;
    fPassedCutsMuon[kRAbsCut]=kFALSE;
    fPassedCutsMuon[kStandardMuonCuts]=kFALSE;

    muonBuffer=0x0;

  }

  Int_t entries=treeMu->GetEntries()-treeEntries;
  cout<<" -> "<<entries<<" muons correctly retrieved."<<endl;

  PostData(1,fOutput);

  return;
}

void AliAnalysisTaskUpsilonTreeDump::Terminate(Option_t *) {
  fOutput=dynamic_cast<TList*>(GetOutputData(1));

  TTree* treeMu=((TTree*)fOutput->At(0));

  cout << "**********************" << endl;
  cout << "* Analysis completed *" << endl;
  cout << "**********************" << endl;
  
  return;
}
