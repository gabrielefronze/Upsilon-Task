#include <iostream>
using std::cout;
using std::endl;

#include "AliMCEvent.h"
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

AliAnalysisTaskRatiosSparse::AliAnalysisTaskRatiosSparse(Bool_t isMC) :
  AliAnalysisTaskSE(""),
  fOutput(0x0)

{
  fCuts=new AliMuonTrackCuts();
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;
  fIsMC=isMC;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskRatiosSparse::AliAnalysisTaskRatiosSparse(const char *name, AliMuonTrackCuts *cuts, Bool_t isMC) :
  AliAnalysisTaskSE(name),
  fOutput(0x0)
{
  fCuts=cuts;
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;
  fIsMC=isMC;

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

  fCuts->SetIsMC(fIsMC);

  PostData(1,fOutput);
}


void AliAnalysisTaskRatiosSparse::UserExec(Option_t *)
{

  if( !fIsMC ){
  	if( !(InputEvent()->GetFiredTriggerClasses()).Contains("CINT7-B-NOPF-MUFAST") ){
  		//cout<<"-> Rejected because of trigger class selection."<<endl;
  		return;
  	} else {
  		//cout<<"-> Accepted after applying trigger class selection"<<endl;
  	}
  }

  ////cout<<"Loading THnSparse"<<endl;
  THnSparseF *sparse=(THnSparseF*)fOutput->At(0);

	////cout<<"Run:"<<InputEvent()->GetRunNumber()<<" Event:"<<fNEvents++;

	//AliInputEventHandler* eventInputHandler=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()));

  Double_t sparseData[kSparseDimension]={0.,0.,0.,0.,0.};

  // reading how much tracks are stored in the input event and the event centrality
  Int_t ntracks=AliAnalysisMuonUtility::GetNTracks(InputEvent());
  Double_t eventCentrality=0.;

  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  if ( multSelection ) eventCentrality=multSelection->GetMultiplicityPercentile("V0M");

  // loop over the tracks to store only muon ones to obtain every possible dimuon
  AliAODTrack* muonBufferData=0x0;
  AliVParticle* muonBufferMC=0x0;
  AliVParticle* motherBufferMC=0x0;
  for(Int_t itrack=0;itrack<ntracks;itrack++){
    muonBufferData=(AliAODTrack*)AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());

    // // check if the track is seen in the muon tracker
    // if ( ! AliAnalysisMuonUtility::IsMuonTrack(muonBufferData) ) continue;
    //
    // // Rapidity acceptance cut.
    // if ( muonBufferData->Eta()<-4.0 || muonBufferData->Eta()>-2.5) {
    // 	//cout<<"Rejected (Eta)"<<endl;
    // 	continue;
    // }
    //
    // // Rabs cut.
    // if ( muonBufferData->GetRAtAbsorberEnd()<17.6 || muonBufferData->GetRAtAbsorberEnd()>89.5) {
    // 	//cout<<"Rejected (Rabs)"<<endl;
    // 	continue;
    // }

    // is the track  selected via standard muon cuts?
    if ( (!fCuts->IsSelected(muonBufferData))) continue;

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
          //else //cout<<"And its mother is an Upsilon! Double cheers!!!"<<endl;
        }
      }
    }

    //Double_t pt=muonBufferData->Pt();

    sparseData[kRapidity]=(Double_t)-muonBufferData->Y();
    sparseData[kPt]=(Double_t)muonBufferData->Pt();
    sparseData[kCentrality]=(Double_t)eventCentrality;
    //sparseData[kTriggerFlag]=(Double_t)AliAnalysisMuonUtility::GetMatchTrigger(muonBufferData);
    sparseData[kLocalBoard]=(Double_t)AliAnalysisMuonUtility::GetLoCircuit(muonBufferData);

    for(Int_t i=0; i<AliAnalysisMuonUtility::GetMatchTrigger(muonBufferData); i++){
      sparseData[kTriggerFlag]=(Double_t)i+1.;
      sparse->Fill(sparseData);
    }

    muonBufferData=0x0;

  }

  PostData(1,fOutput);

  ////cout<<"########## ANALYSIS DONE! ##########"<<endl;

  return;
}

void AliAnalysisTaskRatiosSparse::Terminate(Option_t *) {

  cout<<"########## TERMINATE! ##########"<<endl;

  TH1::SetDefaultSumw2(kTRUE);

  fOutput=dynamic_cast<TList*>(GetOutputData(1));
  THnSparseF *sparse=(THnSparseF*)fOutput->At(0);

  // number of required rapidity bins
  const Int_t nBins=10;

  // how much rebin is needed?
  const Int_t rebinFactor=8;

  Double_t lowRap=2.5;
  Double_t hiRap=4.;
  Double_t binWidth=(hiRap-lowRap)/nBins;

  TH1D *histosAptBuffer;
  TH1D *histosRatio[nBins];

  for (Int_t i = 0; i < nBins; i++) {
    Double_t binLow=lowRap+i*binWidth;
    sparse->GetAxis(kRapidity)->SetRangeUser(binLow, binLow+binWidth);

    sparse->GetAxis(kTriggerFlag)->SetRangeUser(1., 1.);
    histosAptBuffer=sparse->Projection(kPt,"E"); // this is the Apt distribution
    histosAptBuffer->Sumw2();
    histosAptBuffer->Rebin(rebinFactor);

    sparse->GetAxis(kTriggerFlag)->SetRangeUser(2., 2.);
    histosRatio[i]=sparse->Projection(kPt,"E"); // this is the Lpt distribution but will be divided by Apt
    histosRatio[i]->Sumw2();
    histosRatio[i]->Rebin(rebinFactor);

    (histosRatio[i])->SetName(Form("Histo_ratio_%f-%f",binLow,binLow+binWidth));
    (histosRatio[i])->SetTitle(Form("Histo ratio %f-%f",binLow,binLow+binWidth));
    (histosRatio[i])->Divide(histosAptBuffer); // now this histogram contains the ratio of the two

    delete histosAptBuffer;
  }

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

  TFile *file=new TFile(Form("ResponseFunctions%s.root",((fIsMC)?"_MC":"_data")), "RECREATE");

  for (Int_t i = 0; i < nBins; i++) {
    canv->cd(i+1);
    (histosRatio[i])->Draw();
    file->cd();
    (histosRatio[i])->Write();
  }

  cout << "**********************" << endl;
  cout << "* Analysis completed *" << endl;
  cout << "**********************" << endl;

  return;
}
