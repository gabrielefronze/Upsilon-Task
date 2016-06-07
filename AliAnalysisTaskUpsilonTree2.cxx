#include <iostream>
using std::cout;
using std::endl;

#include "TH1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TTree.h"
#include "THnSparse.h"
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
#include "AliAnalysisTaskUpsilonTree2.h"
#include "AliMuonTrackCuts.h"
#include "AliMCEvent.h"


ClassImp(AliAnalysisTaskUpsilonTree2)

AliAnalysisTaskUpsilonTree2::AliAnalysisTaskUpsilonTree2() :
  AliAnalysisTaskSE(""),
  fOutput(0x0),
  fCuts(0x0),
  fNEvents(0),
  fFillTree(kFALSE),
  fIsMC(kFALSE)
{
  fIsMC=kFALSE;
}


AliAnalysisTaskUpsilonTree2::AliAnalysisTaskUpsilonTree2(const char *name, AliMuonTrackCuts *cuts, Bool_t isMC, Bool_t fillTree) :
  AliAnalysisTaskSE(name),
  fOutput(0x0)
{
  fCuts=cuts;
  fCuts->SetAllowDefaultParams(kTRUE);
  fNEvents=0;
  fIsMC=isMC;
  fFillTree=fillTree;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}


AliAnalysisTaskUpsilonTree2::~AliAnalysisTaskUpsilonTree2()
{
  Info("~AliAnalysisTaskUpsilonTree2","Calling Destructor");
  if(fCuts)fCuts=0x0;
  if(fOutput)delete fOutput;
}


void AliAnalysisTaskUpsilonTree2::NotifyRun()
{
  printf("Setting run number for cuts\n");
  fCuts->SetRun(fInputHandler);
}

void AliAnalysisTaskUpsilonTree2::UserCreateOutputObjects()
{
  fOutput = new TList();
  fOutput->SetOwner();

  Float_t *treeData[kSparseDimension];

  if ( fFillTree ){
    TTree *treeDimu = new TTree("MuonData","TTree containing analysis results");
    treeDimu->Branch("dimuon_momentum",&treeData[kMomentum]);
    treeDimu->Branch("dimuon_transverse_momentum",&treeData[kTransverse]);
    treeDimu->Branch("dimuon_rapidity",&treeData[kRapidity]);
    treeDimu->Branch("dimuon_mass",&treeData[kMass]);
    treeDimu->Branch("highest_muon_transverse_momentum",&treeData[kHighMomentumMuonpt]);
    treeDimu->Branch("lowest_muon_transverse_momentum",&treeData[kLowMomentumMuonpt]);
    treeDimu->Branch("enevt_centrality",&treeData[kCentrality]);

    fOutput->AddAt(treeDimu,0);

  } else {
      Int_t nBins[kSparseDimension];
      nBins[kMomentum]=50;
      nBins[kTransverse]=30;
      nBins[kRapidity]=15;
      nBins[kCentrality]=10;
      nBins[kLowMomentumMuonpt]=30;
      nBins[kHighMomentumMuonpt]=30;
      nBins[kMass]=300;

      Double_t xMins[kSparseDimension];
      xMins[kMomentum]=0.;
      xMins[kTransverse]=0.;
      xMins[kRapidity]=2.5;
      xMins[kCentrality]=0.;
      xMins[kLowMomentumMuonpt]=0.;
      xMins[kHighMomentumMuonpt]=0.;
      xMins[kMass]=0.;

      Double_t xMaxs[kSparseDimension];
      xMaxs[kMomentum]=25.;
      xMaxs[kTransverse]=15.;
      xMaxs[kRapidity]=4.;
      xMaxs[kCentrality]=100.;
      xMaxs[kLowMomentumMuonpt]=15.;
      xMaxs[kHighMomentumMuonpt]=15.;
      xMaxs[kMass]=15.;

      THnSparseF *sparseDimu = new THnSparseF("MuonData", "THnSparse containing analysis results", kSparseDimension, nBins, xMins, xMaxs);
      sparseDimu->GetAxis(kMomentum)->SetTitle("dimuon_momentum");
      sparseDimu->GetAxis(kTransverse)->SetTitle("dimuon_transverse_momentum");
      sparseDimu->GetAxis(kRapidity)->SetTitle("dimuon_rapidity");
      sparseDimu->GetAxis(kCentrality)->SetTitle("dimuon_mass");
      sparseDimu->GetAxis(kLowMomentumMuonpt)->SetTitle("highest_muon_transverse_momentum");
      sparseDimu->GetAxis(kHighMomentumMuonpt)->SetTitle("lowest_muon_transverse_momentum");
      sparseDimu->GetAxis(kMass)->SetTitle("enevt_centrality");

      fOutput->AddAt(sparseDimu,0);
  }

  TH1I *histoNEvents = new TH1I("NEvt","NEvt",1,0.,1.);

  fOutput->AddAt(histoNEvents,1);

  PostData(1,fOutput);
}


void AliAnalysisTaskUpsilonTree2::UserExec(Option_t *)
{

  TH1I *histoNEvents =(TH1I*)fOutput->At(1);
  histoNEvents->Fill(0.5);
  histoNEvents = 0x0;
  TObject *dataContainer = fOutput->At(0);

  TTree * treeDimu=0x0;
  Float_t *treeData;

  THnSparseF *sparseDimu=0x0;
  Double_t *sparseData;

  if ( fFillTree ){
    treeDimu = static_cast<TTree*>(dataContainer);
    treeData = new Float_t[kSparseDimension];
    treeDimu->GetBranch("dimuon_momentum")->SetAddress(&treeData[kMomentum]);
    treeDimu->GetBranch("dimuon_transverse_momentum")->SetAddress(&treeData[kTransverse]);
    treeDimu->GetBranch("dimuon_rapidity")->SetAddress(&treeData[kRapidity]);
    treeDimu->GetBranch("dimuon_mass")->SetAddress(&treeData[kMass]);
    treeDimu->GetBranch("highest_muon_transverse_momentum")->SetAddress(&treeData[kHighMomentumMuonpt]);
    treeDimu->GetBranch("lowest_muon_transverse_momentum")->SetAddress(&treeData[kLowMomentumMuonpt]);
    treeDimu->GetBranch("enevt_centrality")->SetAddress(&treeData[kCentrality]);
  } else {
    sparseData = new Double_t[kSparseDimension];
    sparseDimu = static_cast<THnSparseF*>(dataContainer);
  }

	cout<<"Run:"<<InputEvent()->GetRunNumber()<<" Event:"<<fNEvents++;

  // reading how much tracks are stored in the input event and the event centrality
  Int_t ntracks=AliAnalysisMuonUtility::GetNTracks(InputEvent());
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection"));
  Float_t eventCentrality=0.;
  if ( multSelection ) eventCentrality=multSelection->GetMultiplicityPercentile("V0M");

  // array containing all the found muons (eventually after cut apply)
  TObjArray *muonArray = new TObjArray();

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
      //cout<<"Filling data"<<endl;
      if ( fFillTree ) {
        treeData[kMomentum]=dimuon.P();
        treeData[kTransverse]=dimuon.Pt();
        treeData[kRapidity]=dimuon.Rapidity();
        treeData[kMass]=dimuon.M();
        treeData[kHighMomentumMuonpt]=highestMuonp;
        treeData[kLowMomentumMuonpt]=lowestMuonp;
        treeData[kCentrality]=(Double_t)eventCentrality;
        if ( treeData[kMass]>2.5 ) treeDimu->Fill();
      } else {
        sparseData[kMomentum]=dimuon.P();
        sparseData[kTransverse]=dimuon.Pt();
        sparseData[kRapidity]=dimuon.Rapidity();
        sparseData[kMass]=dimuon.M();
        sparseData[kHighMomentumMuonpt]=highestMuonp;
        sparseData[kLowMomentumMuonpt]=lowestMuonp;
        sparseData[kCentrality]=(Double_t)eventCentrality;
        if ( sparseData[kMass]>2.5 )sparseDimu->Fill(sparseData);
      }
      secondMuon=0x0;
    }

    firstMuon=0x0;
  }
  delete muonArray;
  if ( sparseData ) delete[] sparseData;
  if ( treeData ) delete[] treeData;
  dataContainer=0x0;
  treeDimu=0x0;
  sparseDimu=0x0;

  //cout<<"-> "<<treeDimu->GetEntries()-treeEntries<<" dimuons generated and saved."<<endl;

  PostData(1,fOutput);

  return;
}

void AliAnalysisTaskUpsilonTree2::Terminate(Option_t *) {
  fOutput=static_cast<TList*>(GetOutputData(1));

  TObject *dataContainer = fOutput->At(0);
  TTree * treeDimu=0x0;
  THnSparseF *sparseDimu=0x0;

  TCanvas *canvHistos = new TCanvas("canvHistos","canvHistos");
  canvHistos->Divide(2,3);

  //################################################
  // TTREE LIKE
  //################################################
  if ( fFillTree ){
    treeDimu = static_cast<TTree*>(dataContainer);

  ///////////////////////////////////////////////////////////////////
    canvHistos->cd(1)->SetLogy();
    TH1D *histoInvariantMass090FullRap = new TH1D("0%-90% 2.5<eta<4.0","0%-90% 2.5<eta<4.0",600,0.,15.);
   	treeDimu->Project("0%-90% 2.5<eta<4.0","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5");
    histoInvariantMass090FullRap->ShowPeaks();
    histoInvariantMass090FullRap->DrawClone();

  ///////////////////////////////////////////////////////////////////
    canvHistos->cd(2);
    TH1D *histoInvariantMass = new TH1D("All entries","All entries",600,0.,15.);
    treeDimu->Project("All entries","dimuon_mass");
    histoInvariantMass->ShowPeaks();
    histoInvariantMass->DrawClone();

  ///////////////////////////////////////////////////////////////////
  	canvHistos->cd(3)->SetLogy();
    TH1D *histoInvariantMass020FullRap = new TH1D("0%-20% 2.5<eta<4.0","0%-20% 2.5<eta<4.0",600,0.,15.);
   	treeDimu->Project("0%-20% 2.5<eta<4.0","dimuon_mass","enevt_centrality>0 && enevt_centrality<20 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5");
    histoInvariantMass020FullRap->ShowPeaks();
    histoInvariantMass020FullRap->DrawClone();

  ///////////////////////////////////////////////////////////////////
  	canvHistos->cd(4)->SetLogy();
    TH1D *histoInvariantMass2090FullRap = new TH1D("20%-90% 2.5<eta<4.0","20%-90% 2.5<eta<4.0",600,0.,15.);
   	treeDimu->Project("20%-90% 2.5<eta<4.0","dimuon_mass","enevt_centrality>20 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5");
    histoInvariantMass2090FullRap->ShowPeaks();
    histoInvariantMass2090FullRap->DrawClone();

  ///////////////////////////////////////////////////////////////////
  	canvHistos->cd(5)->SetLogy();
    TH1D *histoInvariantMass090SecondRap = new TH1D("0%-90% 2.5<eta<3.2","0%-90% 2.5<eta<3.2",600,0.,15.);
   	treeDimu->Project("0%-90% 2.5<eta<3.2","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-3.2 && dimuon_rapidity<-2.5");
    histoInvariantMass090SecondRap->ShowPeaks();
    histoInvariantMass090SecondRap->DrawClone();

  ///////////////////////////////////////////////////////////////////
  	canvHistos->cd(6)->SetLogy();
    TH1D *histoInvariantMass090FirstRap = new TH1D("0%-90% 3.2<eta<4.0","0%-90% 3.2<eta<4.0",600,0.,15.);
   	treeDimu->Project("0%-90% 3.2<eta<4.0","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-3.2");
    histoInvariantMass090FirstRap->ShowPeaks();
    histoInvariantMass090FirstRap->DrawClone();

  //################################################
  // SPARSE LIKE
  //################################################
  } else {
    sparseDimu = static_cast<THnSparseF*>(dataContainer);

    ///////////////////////////////////////////////////////////////////
      canvHistos->cd(1)->SetLogy();
      TH1D *histoInvariantMass090FullRap = new TH1D();
      histoInvariantMass090FullRap->SetTitle("0%-90% 2.5<eta<4.0");
      sparseDimu->GetAxis(kCentrality)->SetRangeUser(0.,90.);
      histoInvariantMass090FullRap = sparseDimu->Projection(kMass,"E");
      treeDimu->Project("0%-90% 2.5<eta<4.0","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5");
      histoInvariantMass090FullRap->ShowPeaks();
      histoInvariantMass090FullRap->DrawClone();

    ///////////////////////////////////////////////////////////////////
      canvHistos->cd(2);
      TH1D *histoInvariantMass = new TH1D();
      histoInvariantMass->SetTitle("All entries");
      sparseDimu->GetAxis(kCentrality)->SetRangeUser(0.,100.);
      histoInvariantMass = sparseDimu->Projection(kMass,"E");
      treeDimu->Project("All entries","dimuon_mass");
      histoInvariantMass->ShowPeaks();
      histoInvariantMass->DrawClone();

    ///////////////////////////////////////////////////////////////////
      canvHistos->cd(3)->SetLogy();
      TH1D *histoInvariantMass020FullRap = new TH1D();
      histoInvariantMass020FullRap->SetTitle("0%-20% 2.5<eta<4.0");
      sparseDimu->GetAxis(kCentrality)->SetRangeUser(0.,20.);
      histoInvariantMass020FullRap = sparseDimu->Projection(kMass,"E");
      histoInvariantMass020FullRap->ShowPeaks();
      histoInvariantMass020FullRap->DrawClone();

    ///////////////////////////////////////////////////////////////////
      canvHistos->cd(4)->SetLogy();
      TH1D *histoInvariantMass2090FullRap = new TH1D();
      histoInvariantMass2090FullRap->SetTitle("20%-90% 2.5<eta<4.0");
      sparseDimu->GetAxis(kCentrality)->SetRangeUser(20.,90.);
      histoInvariantMass2090FullRap = sparseDimu->Projection(kMass,"E");
      histoInvariantMass2090FullRap->ShowPeaks();
      histoInvariantMass2090FullRap->DrawClone();

    ///////////////////////////////////////////////////////////////////
      canvHistos->cd(5)->SetLogy();
      TH1D *histoInvariantMass090SecondRap = new TH1D();
      histoInvariantMass090SecondRap->SetTitle("0%-90% 2.5<eta<3.2");
      sparseDimu->GetAxis(kCentrality)->SetRangeUser(0.,90.);
      sparseDimu->GetAxis(kRapidity)->SetRangeUser(2.5,3.2);
      histoInvariantMass090SecondRap = sparseDimu->Projection(kMass,"E");
      histoInvariantMass090SecondRap->ShowPeaks();
      histoInvariantMass090SecondRap->DrawClone();

    ///////////////////////////////////////////////////////////////////
      canvHistos->cd(6)->SetLogy();
      TH1D *histoInvariantMass090FirstRap = new TH1D();
      histoInvariantMass090FirstRap->SetTitle("0%-90% 3.2<eta<4.0");
      sparseDimu->GetAxis(kCentrality)->SetRangeUser(0.,90.);
      sparseDimu->GetAxis(kRapidity)->SetRangeUser(3.2,4.0);
      histoInvariantMass090FirstRap = sparseDimu->Projection(kMass,"E");
      histoInvariantMass090FirstRap->ShowPeaks();
      histoInvariantMass090FirstRap->DrawClone();
  }

  dataContainer=0x0;
  treeDimu=0x0;
  sparseDimu=0x0;

  cout << "**********************" << endl;
  cout << "* Analysis completed *" << endl;
  cout << "**********************" << endl;

  return;
}
