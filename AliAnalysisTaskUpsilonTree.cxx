#include <iostream>
using std::cout;
using std::endl;

#include "TH1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TTree.h"
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
#include "AliAnalysisTaskUpsilonTree.h"
#include "AliMuonTrackCuts.h"
#include "AliMCEvent.h"


ClassImp(AliAnalysisTaskUpsilonTree
)

AliAnalysisTaskUpsilonTree::AliAnalysisTaskUpsilonTree() :
  AliAnalysisTaskSE(""),
  fOutput(0x0),
  fCuts(0x0),
  fNEvents(0),
  fTreeData(0x0)
{
  fIsMC=kFALSE;
}


AliAnalysisTaskUpsilonTree::AliAnalysisTaskUpsilonTree(const char *name, AliMuonTrackCuts *cuts, Bool_t isMC) :
  AliAnalysisTaskSE(name),
  fOutput(0x0)
{
  fCuts=cuts;
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;
  fTreeData=(Float_t*)malloc(kSparseDimension*sizeof(Float_t));
  fIsMC=isMC;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskUpsilonTree::~AliAnalysisTaskUpsilonTree()
{
  Info("~AliAnalysisTaskUpsilonTree","Calling Destructor");
  if(fCuts)fCuts=0x0;
  if(fOutput)delete fOutput;
  if(fTreeData)delete fTreeData;
}


void AliAnalysisTaskUpsilonTree::NotifyRun()
{
  printf("Setting run number for cuts\n");
  fCuts->SetRun(fInputHandler);
}

void AliAnalysisTaskUpsilonTree::UserCreateOutputObjects()
{
  fOutput = new TList();
  fOutput->SetOwner();

  TTree *treeDimu=new TTree("MuonData","TTree containing analysis results");
  treeDimu->Branch("dimuon_momentum",&fTreeData[kMomentum]);
  treeDimu->Branch("dimuon_transverse_momentum",&fTreeData[kTransverse]);
  treeDimu->Branch("dimuon_rapidity",&fTreeData[kRapidity]);
  treeDimu->Branch("dimuon_mass",&fTreeData[kMass]);
  treeDimu->Branch("highest_muon_transverse_momentum",&fTreeData[kHighMomentumMuonpt]);
  treeDimu->Branch("lowest_muon_transverse_momentum",&fTreeData[kLowMomentumMuonpt]);
  treeDimu->Branch("enevt_centrality",&fTreeData[kCentrality]);

  TH1I *histoNEvents = new TH1I("NEvt","NEvt",1,0.,1.);

  fOutput->AddAt(treeDimu,0);
  fOutput->AddAt(histoNEvents,1);

  PostData(1,fOutput);
}


void AliAnalysisTaskUpsilonTree::UserExec(Option_t *)
{

  TH1I *histoNEvents =(TH1I*)fOutput->At(1);
  histoNEvents->Fill(0.5);
  histoNEvents = 0x0;
  // sparse that contains data about the dimuon data
  TTree* treeDimu=((TTree*)fOutput->At(0));
  treeDimu->GetBranch("dimuon_momentum")->SetAddress(&fTreeData[kMomentum]);
  treeDimu->GetBranch("dimuon_transverse_momentum")->SetAddress(&fTreeData[kTransverse]);
  treeDimu->GetBranch("dimuon_rapidity")->SetAddress(&fTreeData[kRapidity]);
  treeDimu->GetBranch("dimuon_mass")->SetAddress(&fTreeData[kMass]);
  treeDimu->GetBranch("highest_muon_transverse_momentum")->SetAddress(&fTreeData[kHighMomentumMuonpt]);
  treeDimu->GetBranch("lowest_muon_transverse_momentum")->SetAddress(&fTreeData[kLowMomentumMuonpt]);
  treeDimu->GetBranch("enevt_centrality")->SetAddress(&fTreeData[kCentrality]);

	cout<<"Run:"<<InputEvent()->GetRunNumber()<<" Event:"<<fNEvents++;

  fTreeData=(Float_t*)malloc(kSparseDimension*sizeof(Float_t));

	AliInputEventHandler* eventInputHandler=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()));
	//cout<<eventInputHandler->IsEventSelected()<<endl;
/*	Bool_t isSelected = (eventInputHandler->IsEventSelected();// && AliVEvent::kMUON);
	if ( !isSelected ){
		cout<<"-> Rejected because of physics selection."<<endl;
		return;
	} else cout<<"-> Accepted after applying physics selection";*/

	// if( !(InputEvent()->GetFiredTriggerClasses()).Contains("CMUL7-B-NOPF-MUFAST") ){
	// 	cout<<"-> Rejected because of trigger class selection."<<endl;
	// 	return;
	// } else {
	// 	cout<<"-> Accepted after applying trigger class selection";
	// }

  // reading how much tracks are stored in the input event and the event centrality
  Int_t ntracks=AliAnalysisMuonUtility::GetNTracks(InputEvent());
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  Float_t eventCentrality=0.;
  if ( multSelection ) eventCentrality=multSelection->GetMultiplicityPercentile("V0M");

  // array containing all the found muons (eventually after cut apply)
  TObjArray *muonArray=new TObjArray();

  // loop over the tracks to store only muon ones to obtain every possible dimuon
  AliAODTrack* muonBuffer=0x0;
  AliVParticle* muonBufferMC=0x0;
  AliVParticle* motherBufferMC=0x0;
  for(Int_t itrack=0;itrack<ntracks;itrack++){
    muonBuffer=(AliAODTrack*)AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());

    if ( ! fCuts->IsSelected(muonBuffer) ) continue;

    if( fIsMC ){ // if the analysed run is a MC run the macro excludes any particle not recognized as a muon from upsilon
      Int_t mcDaughterIndex=muonBuffer->GetLabel();
      Int_t mcMotherIndex=0;

      // check if there's any MC truth info related to the particle iTrack
      if ( mcDaughterIndex<0 ) continue; // this particle has no MC truth information
      else {
        muonBufferMC=MCEvent()->GetTrack(mcDaughterIndex);

        // check for particle identity
        if ( TMath::Abs(muonBufferMC->PdgCode())!=13 ) continue; // this particle is not a muon
        else {
          cout<<"It's a muon! Cheers!!!"<<endl;
          mcMotherIndex=muonBufferMC->GetMother();
          if ( mcMotherIndex<0 ) continue;
          motherBufferMC=MCEvent()->GetTrack(mcMotherIndex);

          // check for particle mother identity
          if ( motherBufferMC->PdgCode()!=553 ) continue; // the mother of the studied particle is not a upsilon
          else cout<<"And its mother is an Upsilon! Double cheers!!!"<<endl;
        }
      }
    }

    muonArray->Add((TObject*)muonBuffer);

    muonBuffer=0x0;
  }

  Int_t entries=muonArray->GetEntries();
  Int_t treeEntries=treeDimu->GetEntries();

  cout<<"##### "<<entries<<" muons correctly retrieved";

  AliAODTrack *firstMuon=0x0;
  AliAODTrack *secondMuon=0x0;
  for(Int_t iFirstMuon=0; iFirstMuon<entries-1; iFirstMuon++){
    firstMuon=(AliAODTrack*)muonArray->At(iFirstMuon);
    cout<<"first muon "<<iFirstMuon<<endl;
    Double_t highestMuonp=0.;
    Double_t lowestMuonp=0.;

    highestMuonp=firstMuon->Pt();

    for(Int_t iSecondMuon=iFirstMuon+1; iSecondMuon<entries; iSecondMuon++){
      secondMuon=(AliAODTrack*)muonArray->At(iSecondMuon);

      if (firstMuon->Charge()==secondMuon->Charge()) {
      	//cout<<"Rejected (like-sign)"<<endl;
      	continue;
      }

      //cout<<"second muon "<<iSecondMuon<<endl;
      if(secondMuon->Pt()>highestMuonp){
        lowestMuonp=highestMuonp;
        highestMuonp=secondMuon->Pt();
      } else {
        lowestMuonp=secondMuon->Pt();
      }

      TLorentzVector dimuon=AliAnalysisMuonUtility::GetTrackPair(firstMuon,secondMuon);
      fTreeData[kMomentum]=dimuon.P();
      fTreeData[kTransverse]=dimuon.Pt();
      fTreeData[kRapidity]=dimuon.Rapidity();
      fTreeData[kMass]=dimuon.M();
      fTreeData[kHighMomentumMuonpt]=highestMuonp;
      fTreeData[kLowMomentumMuonpt]=lowestMuonp;
		  fTreeData[kCentrality]=(Double_t)eventCentrality;

      cout<<"FUUUUUUU"<<endl;

      //cout<<"Filling data"<<endl;
      if(fTreeData[kMass]>2.5)treeDimu->Fill();
      cout<<"filled!"<<endl;

     	fTreeData[kMomentum]=0.;
      fTreeData[kTransverse]=0.;
      fTreeData[kRapidity]=0.;
      fTreeData[kMass]=0.;
      fTreeData[kHighMomentumMuonpt]=0.;
      fTreeData[kLowMomentumMuonpt]=0.;
		  fTreeData[kCentrality]=0.;
      secondMuon=0x0;
    }

    firstMuon=0x0;
  }
  delete muonArray;

  //cout<<"-> "<<treeDimu->GetEntries()-treeEntries<<" dimuons generated and saved."<<endl;

  PostData(1,fOutput);

  return;
}

void AliAnalysisTaskUpsilonTree::Terminate(Option_t *) {
  fOutput=dynamic_cast<TList*>(GetOutputData(1));

  TTree* treeDimu=((TTree*)fOutput->At(0));
/*  treeDimu->GetBranch("dimuon_momentum")->SetAddress(&fTreeData[kMomentum]);
  treeDimu->GetBranch("dimuon_transverse_momentum")->SetAddress(&fTreeData[kTransverse]);
  treeDimu->GetBranch("dimuon_rapidity")->SetAddress(&fTreeData[kRapidity]);
  treeDimu->GetBranch("dimuon_mass")->SetAddress(&fTreeData[kMass]);
  treeDimu->GetBranch("highest_muon_transverse_momentum")->SetAddress(&fTreeData[kHighMomentumMuonpt]);
  treeDimu->GetBranch("lowest_muon_transverse_momentum")->SetAddress(&fTreeData[kLowMomentumMuonpt]);
  treeDimu->GetBranch("enevt_centrality")->SetAddress(&fTreeData[kCentrality]);*/

  TCanvas *canvHistos=new TCanvas("canvHistos","canvHistos");
  canvHistos->Divide(2,3);
/*  TH1D *massHisto=new TH1D("histo mass spectrum","histo mass spectrum",1000,2.5,17.);
  treeDimu->Project("histo mass spectrum","dimuon_mass","enevt_centrality>0 & enevt_centrality<10");
  massHisto->Draw();*/

  //TFile *fin=new TFile("RunRoberta.root","READ");
  //TH1D *histoRoberta;
  //fin->GetObject("hMassOS_2m",histoRoberta);
  //histoRoberta->SetLineColor(kRed);
 // histoRoberta->SetLineStyle(3);
  //histoRoberta->Rebin(2);

///////////////////////////////////////////////////////////////////
  canvHistos->cd(1)->SetLogy();
  TH1D *histoInvariantMass090FullRap=new TH1D("0%-90% 2.5<eta<4.0","0%-90% 2.5<eta<4.0",600,0.,15.);
 	treeDimu->Project("0%-90% 2.5<eta<4.0","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5");
  histoInvariantMass090FullRap->ShowPeaks();
  histoInvariantMass090FullRap->DrawClone();
  //histoRoberta->DrawClone("same");

  //fin->Close();

///////////////////////////////////////////////////////////////////
  canvHistos->cd(2);
/*  TGraph *graphRelative=new TGraph();
  graphRelative->SetTitle("Relative variations");
  graphRelative->SetName("Relative variations");
  Double_t binWidth=histoRoberta->GetBinWidth(10);
  Int_t tara=0;
  for(Int_t iBin=1; iBin<histoInvariantMass090FullRap->GetNbinsX(); iBin++){
    Double_t mine=histoInvariantMass090FullRap->GetBinContent(iBin);
    Double_t roberta=histoRoberta->GetBinContent(iBin);
    if (!(iBin*binWidth>2.5)) {
      tara++;
      continue;
    }
    if (roberta!=0.) graphRelative->SetPoint(iBin-1-tara,iBin*binWidth,-(mine-roberta)/roberta*100);
    else graphRelative->SetPoint(iBin-1-tara,iBin*binWidth,0.);

  }
  graphRelative->Draw("ALP");
  graphRelative->SetMarkerStyle(4);
  graphRelative->SetLineColor(kGreen+2);
  graphRelative->SetLineStyle(3);
  graphRelative->SetMarkerColor(kGreen+3);
  graphRelative->SetMarkerSize(0.15);*/
  TH1D *histoInvariantMass=new TH1D("All entries","All entries",600,0.,15.);
  treeDimu->Project("All entries","dimuon_mass");
  histoInvariantMass->ShowPeaks();
  histoInvariantMass->DrawClone();

///////////////////////////////////////////////////////////////////
	canvHistos->cd(3)->SetLogy();
  TH1D *histoInvariantMass020FullRap=new TH1D("0%-20% 2.5<eta<4.0","0%-20% 2.5<eta<4.0",600,0.,15.);
 	treeDimu->Project("0%-20% 2.5<eta<4.0","dimuon_mass","enevt_centrality>0 && enevt_centrality<20 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5");
  histoInvariantMass020FullRap->ShowPeaks();
  histoInvariantMass020FullRap->DrawClone();

///////////////////////////////////////////////////////////////////
	canvHistos->cd(4)->SetLogy();
  TH1D *histoInvariantMass2090FullRap=new TH1D("20%-90% 2.5<eta<4.0","20%-90% 2.5<eta<4.0",600,0.,15.);
 	treeDimu->Project("20%-90% 2.5<eta<4.0","dimuon_mass","enevt_centrality>20 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5");
  histoInvariantMass2090FullRap->ShowPeaks();
  histoInvariantMass2090FullRap->DrawClone();

///////////////////////////////////////////////////////////////////
	canvHistos->cd(5)->SetLogy();
  TH1D *histoInvariantMass090SecondRap=new TH1D("0%-90% 2.5<eta<3.2","0%-90% 2.5<eta<3.2",600,0.,15.);
 	treeDimu->Project("0%-90% 2.5<eta<3.2","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-3.2 && dimuon_rapidity<-2.5");
  histoInvariantMass090SecondRap->ShowPeaks();
  histoInvariantMass090SecondRap->DrawClone();

///////////////////////////////////////////////////////////////////
	canvHistos->cd(6)->SetLogy();
  TH1D *histoInvariantMass090FirstRap=new TH1D("0%-90% 3.2<eta<4.0","0%-90% 3.2<eta<4.0",600,0.,15.);
 	treeDimu->Project("0%-90% 3.2<eta<4.0","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-3.2");
  histoInvariantMass090FirstRap->ShowPeaks();
  histoInvariantMass090FirstRap->DrawClone();

  cout << "**********************" << endl;
  cout << "* Analysis completed *" << endl;
  cout << "**********************" << endl;

  return;
}
