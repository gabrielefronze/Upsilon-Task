#include <iostream>
using std::cout;
using std::endl;

#include "TH1.h"
#include "THnSparse.h"
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
#include "AliAnalysisTaskRatiosSparse.h"
#include "AliMuonTrackCuts.h"


ClassImp(AliAnalysisTaskRatiosSparse)

AliAnalysisTaskRatiosSparse::AliAnalysisTaskRatiosSparse() :
  AliAnalysisTaskSE(""),
  fOutput(0x0)

{
  fCuts=new AliMuonTrackCuts();
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskRatiosSparse::AliAnalysisTaskRatiosSparse(const char *name, AliMuonTrackCuts *cuts) :
  AliAnalysisTaskSE(name),
  fOutput(0x0)
{
  fCuts=cuts;
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskRatiosSparse::~AliAnalysisTaskRatiosSparse()
{
  Info("~AliAnalysisTaskRatiosSparse","Calling Destructor");
  if(fCuts)fCuts=0x0;
  if(fOutput)delete fOutput;
}


void AliAnalysisTaskRatiosSparse::NotifyRun()
{
  printf("Setting run number for cuts\n");
  fCuts->SetRun(fInputHandler);
}

void AliAnalysisTaskRatiosSparse::UserCreateOutputObjects()
{
  fOutput = new TList();
  fOutput->SetOwner();

  const Int_t nBins[kSparseDimension]={   4,  1000, 10,  100,   234};
  const Double_t xMin[kSparseDimension]={0.5,    0., 0.,   0.,   0.5};
  const Double_t xMax[kSparseDimension]={4.5,   12., 4., 100., 234.5};

  THnSparseF *sparse=new THnSparseF("data_sparse","data_sparse",kSparseDimension,nBins,xMin,xMax);
  sparse->GetAxis(kTriggerFlag)->SetTitle("Trigger flag");
  sparse->GetAxis(kPt)->SetTitle("Pt");
  sparse->GetAxis(kRapidity)->SetTitle("Rapidity");
  sparse->GetAxis(kCentrality)->SetTitle("Event centrality");
  sparse->GetAxis(kLocalBoard)->SetTitle("Fired LB");

  fOutput->AddAt(sparse,0);

  PostData(1,fOutput);
}


void AliAnalysisTaskRatiosSparse::UserExec(Option_t *)
{
  THnSparseF *sparse=(THnSparseF*)fOutput->At(0);

	cout<<"Run:"<<InputEvent()->GetRunNumber()<<" Event:"<<fNEvents++;

	AliInputEventHandler* eventInputHandler=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()));

  Double_t sparseData[kSparseDimension]={0.,0.,0.,0.,0.};

	if( !(InputEvent()->GetFiredTriggerClasses()).Contains("CINT7-B-NOPF-MUFAST") ){
		cout<<"-> Rejected because of trigger class selection."<<endl;
		return;
	} else {
		cout<<"-> Accepted after applying trigger class selection"<<endl;
	}

  // reading how much tracks are stored in the input event and the event centrality
  Int_t ntracks=AliAnalysisMuonUtility::GetNTracks(InputEvent());
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  Double_t eventCentrality=multSelection->GetMultiplicityPercentile("V0M");

  // loop over the tracks to store only muon ones to obtain every possible dimuon
  AliAODTrack* muonBuffer=0x0;
  for(Int_t itrack=0;itrack<ntracks;itrack++){
    muonBuffer=(AliAODTrack*)AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());

    // check if the track is seen in the muon tracker
    if ( ! AliAnalysisMuonUtility::IsMuonTrack(muonBuffer) ) continue;

    // Rapidity acceptance cut.
    if ( muonBuffer->Eta()<-4.0 || muonBuffer->Eta()>-2.5) {
    	cout<<"Rejected (Eta)"<<endl;
    	continue;
    }

    // Rabs cut.
    if ( muonBuffer->GetRAtAbsorberEnd()<17.6 || muonBuffer->GetRAtAbsorberEnd()>89.5) {
    	cout<<"Rejected (Rabs)"<<endl;
    	continue;
    }

    // is the track  selected via standard muon cuts?
    if ( ! fCuts->IsSelected(muonBuffer) ) continue;

    Double_t pt=muonBuffer->Pt();

    sparseData[kRapidity]=(Double_t)-muonBuffer->Eta();
    sparseData[kPt]=(Double_t)muonBuffer->Pt();
    sparseData[kCentrality]=(Double_t)eventCentrality;
    //sparseData[kTriggerFlag]=(Double_t)AliAnalysisMuonUtility::GetMatchTrigger(muonBuffer);
    sparseData[kLocalBoard]=(Double_t)AliAnalysisMuonUtility::GetLoCircuit(muonBuffer);

    for(Int_t i=0; i<AliAnalysisMuonUtility::GetMatchTrigger(muonBuffer); i++){
      sparseData[kTriggerFlag]=(Double_t)i+1.;
      sparse->Fill(sparseData);
    }

    muonBuffer=0x0;

  }

  PostData(1,fOutput);

  //cout<<"########## ANALYSIS DONE! ##########"<<endl;

  return;
}

void AliAnalysisTaskRatiosSparse::Terminate(Option_t *) {

  cout<<"########## TERMINATE! ##########"<<endl;

  fOutput=dynamic_cast<TList*>(GetOutputData(1));

  TCanvas *canv=new TCanvas("canv","canv");

  TH1::SetDefaultSumw2(kTRUE);


  cout << "**********************" << endl;
  cout << "* Analysis completed *" << endl;
  cout << "**********************" << endl;

  return;
}
