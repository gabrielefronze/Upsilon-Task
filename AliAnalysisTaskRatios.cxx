#include <iostream>
using std::cout;
using std::endl;

#include "TH1.h"
#include "TGraph.h"
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
#include "AliAnalysisTaskRatios.h"
#include "AliMuonTrackCuts.h"


ClassImp(AliAnalysisTaskRatios)

AliAnalysisTaskRatios::AliAnalysisTaskRatios() :
  AliAnalysisTaskSE(""),
  fOutput(0x0)

{
  fCuts=new AliMuonTrackCuts();
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskRatios::AliAnalysisTaskRatios(const char *name, AliMuonTrackCuts *cuts) :
  AliAnalysisTaskSE(name),
  fOutput(0x0)
{
  fCuts=cuts;
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskRatios::~AliAnalysisTaskRatios()
{
  Info("~AliAnalysisTaskRatios","Calling Destructor");
  if(fCuts)fCuts=0x0;
  if(fOutput)delete fOutput;
}


void AliAnalysisTaskRatios::NotifyRun()
{
  printf("Setting run number for cuts\n");
  fCuts->SetRun(fInputHandler);
}

void AliAnalysisTaskRatios::UserCreateOutputObjects()
{
  fOutput = new TList();
  fOutput->SetOwner();

  TH1D *histoLpt=new TH1D("Lpt","Lpt",1000,0.,20.);
  TH1D *histoApt=new TH1D("Apt","Apt",1000,0.,20.);

  TH1D *histoLpt253=new TH1D("Lpt253","Lpt253",1000,0.,20.);
  TH1D *histoApt253=new TH1D("Apt253","Apt253",1000,0.,20.);

  TH1D *histoLpt335=new TH1D("Lpt335","Lpt335",1000,0.,20.);
  TH1D *histoApt335=new TH1D("Apt335","Apt335",1000,0.,20.);

  TH1D *histoLpt354=new TH1D("Lpt354","Lpt354",1000,0.,20.);
  TH1D *histoApt354=new TH1D("Apt354","Apt354",1000,0.,20.);

  fOutput->AddAt(histoLpt,0);
  fOutput->AddAt(histoApt,1);

  fOutput->AddAt(histoLpt253,2);
  fOutput->AddAt(histoApt253,3);

  fOutput->AddAt(histoLpt335,4);
  fOutput->AddAt(histoApt335,5);

  fOutput->AddAt(histoLpt354,6);
  fOutput->AddAt(histoApt354,7);

  PostData(1,fOutput);
}


void AliAnalysisTaskRatios::UserExec(Option_t *)
{
  TH1D *histoApt=(TH1D*)fOutput->At(0);
  TH1D *histoLpt=(TH1D*)fOutput->At(1);
  TH1D *histoLpt253=(TH1D*)fOutput->At(2);
  TH1D *histoApt253=(TH1D*)fOutput->At(3);
  TH1D *histoLpt335=(TH1D*)fOutput->At(4);
  TH1D *histoApt335=(TH1D*)fOutput->At(5);
  TH1D *histoLpt354=(TH1D*)fOutput->At(6);
  TH1D *histoApt354=(TH1D*)fOutput->At(7);

	cout<<"Run:"<<InputEvent()->GetRunNumber()<<" Event:"<<fNEvents++;

	AliInputEventHandler* eventInputHandler=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()));

	if( !(InputEvent()->GetFiredTriggerClasses()).Contains("CINT7-B-NOPF-MUFAST") ){
		cout<<"-> Rejected because of trigger class selection."<<endl;
		return;
	} else {
		cout<<"-> Accepted after applying trigger class selection"<<endl;
	}

  // reading how much tracks are stored in the input event and the event centrality
  Int_t ntracks=AliAnalysisMuonUtility::GetNTracks(InputEvent());
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  Float_t eventCentrality=multSelection->GetMultiplicityPercentile("V0M");

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

    Double_t rapidity=-muonBuffer->Eta();
    Double_t pt=muonBuffer->Pt();

    // keep only low-pt matching tracks
    if ( AliAnalysisMuonUtility::GetMatchTrigger(muonBuffer)>1 ){ //if Apt
      cout<<"Apt"<<endl;
      histoApt->Fill(pt);
      if ( rapidity>=2.5 && rapidity<3. ) histoApt253->Fill(pt);
      if ( rapidity>=3. && rapidity<3.5 ) histoApt335->Fill(pt);
      if ( rapidity>=3.5 && rapidity<4. ) histoApt354->Fill(pt);
    }

    if ( AliAnalysisMuonUtility::GetMatchTrigger(muonBuffer)>2 ){ //if Lpt
      cout<<"Lpt"<<endl;
      histoLpt->Fill(pt);
      if ( rapidity>=2.5 && rapidity<3. ) histoLpt253->Fill(pt);
      if ( rapidity>=3. && rapidity<3.5 ) histoLpt335->Fill(pt);
      if ( rapidity>=3.5 && rapidity<4. ) histoLpt354->Fill(pt);
    }

    muonBuffer=0x0;

  }

  PostData(1,fOutput);

  //cout<<"########## ANALYSIS DONE! ##########"<<endl;

  return;
}

void AliAnalysisTaskRatios::Terminate(Option_t *) {

  cout<<"########## TERMINATE! ##########"<<endl;

  fOutput=dynamic_cast<TList*>(GetOutputData(1));

  TCanvas *canv=new TCanvas("canv","canv");

  TH1D *histoApt=(TH1D*)fOutput->At(0);
  TH1D *histoLpt=(TH1D*)fOutput->At(1);
  TH1D *histoLpt253=(TH1D*)fOutput->At(2);
  TH1D *histoApt253=(TH1D*)fOutput->At(3);
  TH1D *histoLpt335=(TH1D*)fOutput->At(4);
  TH1D *histoApt335=(TH1D*)fOutput->At(5);
  TH1D *histoLpt354=(TH1D*)fOutput->At(6);
  TH1D *histoApt354=(TH1D*)fOutput->At(7);
  histoLpt->Rebin(8);
  histoApt->Rebin(8);
  histoLpt253->Rebin(8);
  histoApt253->Rebin(8);
  histoLpt335->Rebin(8);
  histoApt335->Rebin(8);
  histoLpt354->Rebin(8);
  histoApt354->Rebin(8);


  TH1D *histoRatio=new TH1D("histoRatio","histoRatio",1000,0.,20.);
  TH1D *histoRatio253=new TH1D("histoRatio253","histoRatio253",1000,0.,20.);
  TH1D *histoRatio335=new TH1D("histoRatio335","histoRatio335",1000,0.,20.);
  TH1D *histoRatio354=new TH1D("histoRatio354","histoRatio354",1000,0.,20.);
  histoRatio->Rebin(8);
  histoRatio253->Rebin(8);
  histoRatio335->Rebin(8);
  histoRatio354->Rebin(8);
  histoRatio->Divide(histoLpt,histoApt);
  histoRatio253->Divide(histoLpt253,histoApt253);
  histoRatio335->Divide(histoLpt335,histoApt335);
  histoRatio354->Divide(histoLpt354,histoApt354);
  histoRatio->SetTitle("integrated ratio");
  histoRatio253->SetTitle("ratio 2.5-3.0");
  histoRatio335->SetTitle("ratio 3.0-3.5");
  histoRatio354->SetTitle("ratio 3.5-4.0");

  canv->Divide(2,2);
  canv->cd(1);
  histoRatio->GetXaxis()->SetRangeUser(0.,20.);
  histoRatio->Draw();
  canv->cd(2);
  histoRatio253->GetXaxis()->SetRangeUser(0.,20.);
  histoRatio253->Draw();
  canv->cd(3);
  histoRatio335->GetXaxis()->SetRangeUser(0.,20.);
  histoRatio335->Draw();
  canv->cd(4);
  histoRatio354->GetXaxis()->SetRangeUser(0.,20.);
  histoRatio354->Draw();

  cout << "**********************" << endl;
  cout << "* Analysis completed *" << endl;
  cout << "**********************" << endl;

  return;
}
