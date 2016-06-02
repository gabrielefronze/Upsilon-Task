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
#include "THnSparse.h"
#include "TObjArray.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliCentrality.h"
#include "AliAODTrack.h"
#include "TLorentzVector.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMultSelection.h"
#include "AliAnalysisTaskUpsilon.h"
#include "AliMuonTrackCuts.h"


ClassImp(AliAnalysisTaskUpsilon
)

AliAnalysisTaskUpsilon::AliAnalysisTaskUpsilon() :
  AliAnalysisTaskSE(""),
  fAODEvent(0x0),
  fESDEvent(0x0),
  fOutput(0x0)
{
  fCuts=new AliMuonTrackCuts();
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskUpsilon::AliAnalysisTaskUpsilon(const char *name, AliMuonTrackCuts *cuts) :
  AliAnalysisTaskSE(name),
  fAODEvent(0x0),
  fESDEvent(0x0),
  fOutput(0x0)
{
  fCuts=cuts;
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskUpsilon::~AliAnalysisTaskUpsilon()
{
  Info("~AliAnalysisTaskUpsilon","Calling Destructor");
  if(fCuts)fCuts=0x0;
  if(fOutput)
  {
    delete fOutput;
    fOutput = 0x0;
  }
}


void AliAnalysisTaskUpsilon::NotifyRun()
{
  printf("Setting run number for cuts\n");
  fCuts->SetRun(fInputHandler);
}

void AliAnalysisTaskUpsilon::UserCreateOutputObjects()
{
  fOutput = new TList();
  fOutput->SetOwner();

  Int_t bins[kSparseDimension]={17,17,10,10,17,17,340};
  Double_t mins[kSparseDimension]={0.0,0.0,-4.0,000.5,0.0,0.0,0.0};
  Double_t maxs[kSparseDimension]={17.,17.,-2.5,100.5,34.,34.,17.};

  THnSparseD *sparseDimu=new THnSparseD("MuonData","MuonData",kSparseDimension,bins,mins,maxs,1);

  sparseDimu->GetAxis(kMomentum)->SetTitle("Momentum");
  sparseDimu->GetAxis(kTransverse)->SetTitle("Transverse momentum");
  sparseDimu->GetAxis(kRapidity)->SetTitle("Rapidity");
  sparseDimu->GetAxis(kCentrality)->SetTitle("Event centrality");
  sparseDimu->GetAxis(kMass)->SetTitle("Mass");
  sparseDimu->GetAxis(kHighMomentumMuonp)->SetTitle("Momentum of leading muon");
  sparseDimu->GetAxis(kLowMomentumMuonp)->SetTitle("Momentum of second muon");  

  fOutput->AddAt(sparseDimu,0);

  PostData(1,fOutput);
}


void AliAnalysisTaskUpsilon::UserExec(Option_t *)
{
	cout<<"Run:"<<InputEvent()->GetRunNumber()<<" Event:"<<fNEvents++;
	//AliInputEventHandler* eventInputHandler=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()));
	//cout<<eventInputHandler->IsEventSelected()<<endl;
/*	Bool_t isSelected = (eventInputHandler->IsEventSelected());
	if ( !isSelected ){
		cout<<"-> Rejected because of physics selection."<<endl;
		return;
	} else cout<<"-> Accepted after applying physics selection";*/

	if( !(InputEvent()->GetFiredTriggerClasses()).Contains("CMUL7-B-NOPF-MUFAST") ){
		cout<<"-> Rejected because of trigger class selection."<<endl;
		return;
	} else {
		cout<<"-> Accepted after applying trigger class selection";
	}
  // sparse that contains data about the dimuon data
  THnSparseD* sparseDimu=((THnSparseD*)fOutput->At(0));

  // reading how much tracks are stored in the input event and the event centrality
  Int_t ntracks=AliAnalysisMuonUtility::GetNTracks(InputEvent());
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  Float_t eventCentrality=multSelection->GetMultiplicityPercentile("V0M");

  // double array that will be used for the sparse filling
  Double_t sparseData[kSparseDimension];
  sparseData[kCentrality]=(Double_t)eventCentrality;

  // array containing all the found muons (eventually after cut apply)
  TObjArray *muonArray=new TObjArray();

  // loop over the tracks to store only muon ones to obtain every possible dimuon
  AliAODTrack* muonBuffer=0x0;
  for(Int_t itrack=0;itrack<ntracks;itrack++){
    muonBuffer=(AliAODTrack*)AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());

    // check if the track is seen in the muon tracker
    if ( ! AliAnalysisMuonUtility::IsMuonTrack(muonBuffer) ) continue;

    // keep only low-pt matching tracks
    if ( AliAnalysisMuonUtility::GetMatchTrigger(muonBuffer)<2 ) continue;

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
    //if ( ! fCuts->IsSelected(muonBuffer) ) continue;

    muonArray->Add((TObject*)muonBuffer);

    muonBuffer=0x0;
  }

  Int_t entries=muonArray->GetEntries();
  Int_t sparseEntries=sparseDimu->GetEntries();
  cout<<"-> "<<entries<<" muons correctly retrieved";

  AliAODTrack *firstMuon=0x0;
  AliAODTrack *secondMuon=0x0;
  TLorentzVector dimuon;
  for(Int_t iFirstMuon=0; iFirstMuon<entries-1; iFirstMuon++){
    firstMuon=(AliAODTrack*)muonArray->At(iFirstMuon);
    //cout<<"first muon "<<iFirstMuon<<endl;
    Double_t highestMuonp=0.;
    Double_t lowestMuonp=0.;

    highestMuonp=firstMuon->Pt();

    for(Int_t iSecondMuon=iFirstMuon+1; iSecondMuon<entries; iSecondMuon++){
      secondMuon=(AliAODTrack*)muonArray->At(iSecondMuon);

      if (firstMuon->Charge()==secondMuon->Charge()) {
      	cout<<"Rejected (like-sign)"<<endl;
      	continue;
      }

      //cout<<"second muon "<<iSecondMuon<<endl;
      if(secondMuon->Pt()>highestMuonp){
        lowestMuonp=highestMuonp;
        highestMuonp=secondMuon->Pt();
      } else {
        lowestMuonp=secondMuon->Pt();
      }

      dimuon=AliAnalysisMuonUtility::GetTrackPair(firstMuon,secondMuon);
      sparseData[kMomentum]=dimuon.P();
      sparseData[kTransverse]=dimuon.Pt();
      sparseData[kRapidity]=dimuon.Rapidity();
      sparseData[kMass]=dimuon.M();
      sparseData[kHighMomentumMuonp]=highestMuonp;
      sparseData[kLowMomentumMuonp]=lowestMuonp;

      //cout<<"Filling data"<<endl;
      sparseDimu->Fill(sparseData);

      secondMuon=0x0;
    }

    firstMuon=0x0;
  }
  delete muonArray;

  cout<<"-> "<<sparseDimu->GetEntries()-sparseEntries<<" dimuons generated and saved."<<endl;

  PostData(1,fOutput);

  return;
}

void AliAnalysisTaskUpsilon::Terminate(Option_t *) {
  fOutput=dynamic_cast<TList*>(GetOutputData(1));

  TCanvas *canvHistos=new TCanvas("canvHistos","canvHistos");
  canvHistos->Divide(2,3);

  THnSparseD* sparseDimu=((THnSparseD*)fOutput->At(0));
  sparseDimu->GetAxis(kHighMomentumMuonp)->SetRangeUser(4., 100000000.);
  sparseDimu->GetAxis(kLowMomentumMuonp)->SetRangeUser(2., 100000000.);
  sparseDimu->GetAxis(kRapidity)->SetRangeUser(0., 12.);
  sparseDimu->GetAxis(kMomentum)->SetRangeUser(4.,100000000.);

/*  TH1D* lowMomentumMuSpectrum=sparseDimu->Projection(kLowMomentumMuonp);
  TH1D* highMomentumMuSpectrum=sparseDimu->Projection(kHighMomentumMuonp);
  canvHistos->cd(2);
  lowMomentumMuSpectrum->Draw();
  highMomentumMuSpectrum->SetLineColor(kRed);
	highMomentumMuSpectrum->Draw("SAME");*/

///////////////////////////////////////////////////////////////////
  sparseDimu->GetAxis(kCentrality)->SetRangeUser(0., 90.);
  sparseDimu->GetAxis(kRapidity)->SetRangeUser(-4.0, -2.5);

  canvHistos->cd(1);
  TH1D *histoInvariantMass090FullRap=sparseDimu->Projection(kMass);
 	histoInvariantMass090FullRap->Rebin(2);
  histoInvariantMass090FullRap->GetXaxis()->SetRangeUser(2., 12.);
  histoInvariantMass090FullRap->SetTitle("0%-90% -4.0<eta<-2.5");
  histoInvariantMass090FullRap->ShowPeaks();
  histoInvariantMass090FullRap->DrawClone();

///////////////////////////////////////////////////////////////////
  sparseDimu->GetAxis(kCentrality)->SetRangeUser(0., 20.);
  sparseDimu->GetAxis(kRapidity)->SetRangeUser(-4.0, -2.5);

	canvHistos->cd(3);
  TH1D *histoInvariantMass020FullRap=sparseDimu->Projection(kMass);
  histoInvariantMass020FullRap->Rebin(2);
  histoInvariantMass020FullRap->GetXaxis()->SetRangeUser(2., 12.);
  histoInvariantMass020FullRap->SetTitle("0%-20% -4.0<eta<-2.5");
  histoInvariantMass020FullRap->ShowPeaks();
  histoInvariantMass020FullRap->DrawClone();

///////////////////////////////////////////////////////////////////
  sparseDimu->GetAxis(kCentrality)->SetRangeUser(20., 90.);
  sparseDimu->GetAxis(kRapidity)->SetRangeUser(-4.0, -2.5);

	canvHistos->cd(4);
  TH1D *histoInvariantMass2090FullRap=sparseDimu->Projection(kMass);
  histoInvariantMass2090FullRap->Rebin(2);
  histoInvariantMass2090FullRap->GetXaxis()->SetRangeUser(2., 12.);
  histoInvariantMass2090FullRap->SetTitle("20%-90% -4.0<eta<-2.5");
  histoInvariantMass2090FullRap->ShowPeaks();
  histoInvariantMass2090FullRap->DrawClone();

///////////////////////////////////////////////////////////////////
  sparseDimu->GetAxis(kCentrality)->SetRangeUser(0., 90.);
  sparseDimu->GetAxis(kRapidity)->SetRangeUser(-4.0, -3.2);

	canvHistos->cd(5);
  TH1D *histoInvariantMass090FirstRap=sparseDimu->Projection(kMass);
  histoInvariantMass090FirstRap->Rebin(2);
  histoInvariantMass090FirstRap->GetXaxis()->SetRangeUser(2., 12.);
  histoInvariantMass090FirstRap->SetTitle("0%-90% -4.0<eta<-3.2");
  histoInvariantMass090FirstRap->ShowPeaks();
  histoInvariantMass090FirstRap->DrawClone();

///////////////////////////////////////////////////////////////////
  sparseDimu->GetAxis(kCentrality)->SetRangeUser(0., 90.);
  sparseDimu->GetAxis(kRapidity)->SetRangeUser(-3.2, -2.5);

	canvHistos->cd(6);
  TH1D *histoInvariantMass090SecondRap=sparseDimu->Projection(kMass);
  histoInvariantMass090SecondRap->Rebin(2);
  histoInvariantMass090SecondRap->GetXaxis()->SetRangeUser(2., 12.);
  histoInvariantMass090SecondRap->SetTitle("0%-90% -3.2<eta<-2.5");
  histoInvariantMass090SecondRap->ShowPeaks();
  histoInvariantMass090SecondRap->DrawClone();

/*  TH1D *histoDNOverDM=sparseDimu->Projection(kMass);
  for(Int_t i=1; i<histoDNOverDM->GetNbinsX(); i++){
  	printf("bin: %d prev: %f new: %f\n",i,histoDNOverDM->GetBinContent(i),histoDNOverDM->GetBinContent(i)*((Double_t)i*histoDNOverDM->GetBinWidth(i)));
  	histoDNOverDM->SetBinContent(i,histoDNOverDM->GetBinContent(i)*((Double_t)i*histoDNOverDM->GetBinWidth(i))*((Double_t)i*histoDNOverDM->GetBinWidth(i)));
  }
  histoDNOverDM->GetXaxis()->SetRangeUser(2.5, 14.);
  histoDNOverDM->DrawClone();  */
  cout << "**********************" << endl;
  cout << "* Analysis completed *" << endl;
  cout << "**********************" << endl;

  return;
}
