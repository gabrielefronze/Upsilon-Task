#include <iostream>
using std::cout;
using std::endl;
using namespace std;

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
#include "AliAnalysisTaskUpsilonTreeAEff.h"
#include "AliMuonTrackCuts.h"
#include "AliMCEvent.h"


ClassImp(AliAnalysisTaskUpsilonTreeAEff
)

AliAnalysisTaskUpsilonTreeAEff::AliAnalysisTaskUpsilonTreeAEff() :
  AliAnalysisTaskSE(""),
  fOutput(0x0),
  fCuts(0x0),
  fNEvents(0),
  fTreeData(0x0)
{
  fIsMC=kFALSE;
}


AliAnalysisTaskUpsilonTreeAEff::AliAnalysisTaskUpsilonTreeAEff(const char *name, AliMuonTrackCuts *cuts, Bool_t isMC) :
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


AliAnalysisTaskUpsilonTreeAEff::~AliAnalysisTaskUpsilonTreeAEff()
{
  Info("~AliAnalysisTaskUpsilonTreeAEff","Calling Destructor");
  if(fCuts)fCuts=0x0;
  if(fOutput)delete fOutput;
  if(fTreeData)delete fTreeData;
}


void AliAnalysisTaskUpsilonTreeAEff::NotifyRun()
{
  printf("Setting run number for cuts\n");
  fCuts->SetRun(fInputHandler);
}

void AliAnalysisTaskUpsilonTreeAEff::UserCreateOutputObjects()
{
  fOutput = new TList();
  fOutput->SetOwner();

  TTree *treeDimu=new TTree("DimuonData","TTree containing analysis results");
  treeDimu->Branch("dimuon_momentum",&fTreeData[kMomentum]);
  treeDimu->Branch("dimuon_transverse_momentum",&fTreeData[kTransverse]);
  treeDimu->Branch("dimuon_rapidity",&fTreeData[kRapidity]);
  treeDimu->Branch("dimuon_mass",&fTreeData[kMass]);
  treeDimu->Branch("highest_muon_transverse_momentum",&fTreeData[kHighMomentumMuonpt]);
  treeDimu->Branch("lowest_muon_transverse_momentum",&fTreeData[kLowMomentumMuonpt]);
  treeDimu->Branch("enevt_centrality",&fTreeData[kCentrality]);

  TTree *treeDimuMC=new TTree("DimuonDataMC","TTree containing MC results");
  treeDimuMC->Branch("dimuon_momentum",&fTreeData[kMomentum]);
  treeDimuMC->Branch("dimuon_transverse_momentum",&fTreeData[kTransverse]);
  treeDimuMC->Branch("dimuon_rapidity",&fTreeData[kRapidity]);
  treeDimuMC->Branch("dimuon_mass",&fTreeData[kMass]);
  treeDimuMC->Branch("highest_muon_transverse_momentum",&fTreeData[kHighMomentumMuonpt]);
  treeDimuMC->Branch("lowest_muon_transverse_momentum",&fTreeData[kLowMomentumMuonpt]);
  treeDimuMC->Branch("enevt_centrality",&fTreeData[kCentrality]);

  TH1I *histoNEvents = new TH1I("NEvt","NEvt",1,0.,2.);

  fOutput->AddAt(treeDimu,0);
  fOutput->AddAt(histoNEvents,1);
  fOutput->AddAt(treeDimuMC,2);

  PostData(1,fOutput);
}


void AliAnalysisTaskUpsilonTreeAEff::UserExec(Option_t *)
{
  TH1I *histoNEvents =(TH1I*)fOutput->At(1);
  histoNEvents->Fill(1.);

  fTreeData=(Float_t*)malloc(kSparseDimension*sizeof(Float_t));

  for (size_t i = 0; i < kSparseDimension; i++) {
    fTreeData[i]=0.;
  }

  // sparse that contains data about the dimuon data
  TTree* treeDimu=((TTree*)fOutput->At(0));
  treeDimu->GetBranch("dimuon_momentum")->SetAddress(&fTreeData[kMomentum]);
  treeDimu->GetBranch("dimuon_transverse_momentum")->SetAddress(&fTreeData[kTransverse]);
  treeDimu->GetBranch("dimuon_rapidity")->SetAddress(&fTreeData[kRapidity]);
  treeDimu->GetBranch("dimuon_mass")->SetAddress(&fTreeData[kMass]);
  treeDimu->GetBranch("highest_muon_transverse_momentum")->SetAddress(&fTreeData[kHighMomentumMuonpt]);
  treeDimu->GetBranch("lowest_muon_transverse_momentum")->SetAddress(&fTreeData[kLowMomentumMuonpt]);
  treeDimu->GetBranch("enevt_centrality")->SetAddress(&fTreeData[kCentrality]);

  cout<<"################################################"<<endl;
	cout<<"Run:"<<InputEvent()->GetRunNumber()<<" Event:"<<fNEvents++<<endl;

	//AliInputEventHandler* eventInputHandler=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()));
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

  TObjArray muPlus,muMinus;
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

        if ( muonBufferMC->Charge()>0. ) muPlus.AddAt(muonBufferData, mcMotherIndex);
        if ( muonBufferMC->Charge()<0. ) {
          muMinus.AddAt(muonBufferData, mcMotherIndex);
          motherIndexes.push_back (mcMotherIndex);
        }

        //cout<<"found upsilon #"<<mcMotherIndex<<endl;
        //else //cout<<"And its mother is an Upsilon! Double cheers!!!"<<endl;
      }
    }
  }

  AliAODTrack *firstMuon;
  AliAODTrack *secondMuon;
  for(std::vector<Int_t>::iterator mothersIterator = motherIndexes.begin(); mothersIterator != motherIndexes.end(); ++mothersIterator) {
    //cout<<"mother #"<<*mothersIterator<<endl;
    firstMuon=(AliAODTrack*)muPlus.At(*mothersIterator);
    secondMuon=(AliAODTrack*)muMinus.At(*mothersIterator);

    if ( !firstMuon || !secondMuon ) continue;

    TLorentzVector dimuon = AliAnalysisMuonUtility::GetTrackPair(firstMuon,secondMuon);

    fTreeData[kMomentum]=dimuon.P();
    fTreeData[kTransverse]=dimuon.Pt();
    fTreeData[kRapidity]=TMath::Abs(dimuon.Y());
    fTreeData[kMass]=dimuon.M();
    fTreeData[kHighMomentumMuonpt]=0.;
    fTreeData[kLowMomentumMuonpt]=0.;
    fTreeData[kCentrality]=(Double_t)eventCentrality;

    treeDimu->Fill();
    cout<<"Filled!"<<endl;

    cout<<endl;
    cout<<fTreeData[kMomentum]<<endl;
    cout<<fTreeData[kTransverse]<<endl;
    cout<<fTreeData[kRapidity]<<endl;
    cout<<fTreeData[kMass]<<endl;
    cout<<fTreeData[kHighMomentumMuonpt]<<endl;
    cout<<fTreeData[kLowMomentumMuonpt]<<endl;
    cout<<fTreeData[kCentrality]<<endl;
    cout<<endl;
    cout<<endl;
  }

  treeDimu = 0x0;
  firstMuon = 0x0;
  secondMuon = 0x0;


  TTree* treeDimuMC=((TTree*)fOutput->At(2));
  treeDimuMC->GetBranch("dimuon_momentum")->SetAddress(&fTreeData[kMomentum]);
  treeDimuMC->GetBranch("dimuon_transverse_momentum")->SetAddress(&fTreeData[kTransverse]);
  treeDimuMC->GetBranch("dimuon_rapidity")->SetAddress(&fTreeData[kRapidity]);
  treeDimuMC->GetBranch("dimuon_mass")->SetAddress(&fTreeData[kMass]);
  treeDimuMC->GetBranch("highest_muon_transverse_momentum")->SetAddress(&fTreeData[kHighMomentumMuonpt]);
  treeDimuMC->GetBranch("lowest_muon_transverse_momentum")->SetAddress(&fTreeData[kLowMomentumMuonpt]);
  treeDimuMC->GetBranch("enevt_centrality")->SetAddress(&fTreeData[kCentrality]);

  Int_t nMCtracks=MCEvent()->GetNumberOfTracks();
  for (size_t i = 0; i < kSparseDimension; i++) {
    fTreeData[i]=0.;
  }

  //cout<<"On MC from now"<<endl;

  AliVParticle *dimuonBuffer = 0x0;
  AliVParticle *muBuffer = 0x0;
  AliVParticle *firstMuonMC = 0x0;
  AliVParticle *secondMuonMC = 0x0;
  for(Int_t itrackMC=0;itrackMC<nMCtracks;itrackMC++){
    dimuonBuffer=MCEvent()->GetTrack(itrackMC);

    //cout<<(AliAnalysisMuonUtility::GetTrackHistory(dimuonBuffer, MCEvent(), kTRUE)).Data()<<endl;
    //cout<<(MCEvent()->GetGenerator(dimuonBuffer->GetGeneratorIndex())).Data()<<endl;
    if ( dimuonBuffer->PdgCode()!=553 ) continue;
    else cout<<"Upsilon!!!"<<endl;

    //cout<<"OK pippo2"<<endl;
    Int_t firstD = dimuonBuffer->GetFirstDaughter();
    Int_t secondD = dimuonBuffer->GetLastDaughter();

    printf("%d %d %d \n",firstD,secondD,nMCtracks);

    if ( firstD < 0 ) continue;

    for ( Int_t iDau=firstD; iDau<secondD; iDau++ ){
      muBuffer=MCEvent()->GetTrack(iDau);
      if ( TMath::Abs(muBuffer->PdgCode()) == 13 ){
        if ( muBuffer->Charge()>0. ) firstMuonMC = MCEvent()->GetTrack(iDau);
        if ( muBuffer->Charge()<0. ) secondMuonMC = MCEvent()->GetTrack(iDau);
      }
    }

    if ( !firstMuonMC || !secondMuonMC ) continue;

    cout<<firstMuonMC->Eta()<<endl;
    if ( firstMuonMC->Eta()<-4.0 || firstMuonMC->Eta()>-2.5 ) continue;
    cout<<"OK1"<<endl;

    cout<<secondMuonMC->Eta()<<endl;
    if ( secondMuonMC->Eta()<-4.0 || secondMuonMC->Eta()>-2.5 ) continue;
    cout<<"OK2"<<endl;

    // if ( !fCuts->IsSelected(firstMuonMC) || !fCuts->IsSelected(secondMuonMC) ) continue;

    fTreeData[kMomentum]=dimuonBuffer->P();
    fTreeData[kTransverse]=dimuonBuffer->Pt();
    fTreeData[kRapidity]=TMath::Abs(dimuonBuffer->Y());
    fTreeData[kMass]=dimuonBuffer->M();
    fTreeData[kHighMomentumMuonpt]=0.;
    fTreeData[kLowMomentumMuonpt]=0.;
    fTreeData[kCentrality]=(Double_t)eventCentrality;

    treeDimuMC->Fill();

    cout<<endl;
    cout<<fTreeData[kMomentum]<<endl;
    cout<<fTreeData[kTransverse]<<endl;
    cout<<fTreeData[kRapidity]<<endl;
    cout<<fTreeData[kMass]<<endl;
    cout<<fTreeData[kHighMomentumMuonpt]<<endl;
    cout<<fTreeData[kLowMomentumMuonpt]<<endl;
    cout<<fTreeData[kCentrality]<<endl;
    cout<<endl;
    cout<<endl;

    cout<<"filled MC!"<<endl;
    firstMuonMC = 0x0;
    secondMuonMC = 0x0;
  }

  //cout<<"-> "<<treeDimu->GetEntries()-treeEntries<<" dimuons generated and saved."<<endl;

  PostData(1,fOutput);

  return;
}

void AliAnalysisTaskUpsilonTreeAEff::Terminate(Option_t *) {
  fOutput=dynamic_cast<TList*>(GetOutputData(1));

  TH1I *histoNEvents =(TH1I*)fOutput->At(1);
  cout<<histoNEvents->IsA()->GetName()<<endl;
  Double_t nEvents = (Double_t)histoNEvents->GetBinContent(1);
  cout<<nEvents<<endl;

  TCanvas *canvHistos=new TCanvas("canvHistos","canvHistos");
  canvHistos->Divide(2,1);

  canvHistos->cd(1);
  TTree* treeDimu=((TTree*)fOutput->At(0));
  cout<<treeDimu->GetEntries()<<endl;
  TH1D *rapiditySpectrum1 = new TH1D("rapidity_distribution_2.5-4.0_Data","Rapidity distribution 2.5<#eta<4.0",3,2.5,4.0);
  // rapiditySpectrum1->Sumw2(kTRUE);
  treeDimu->Project("rapidity_distribution_2.5-4.0_Data","dimuon_rapidity","dimuon_rapidity>2.5 && dimuon_rapidity<4.0");
  rapiditySpectrum1->Draw("E");

  canvHistos->cd(1);
  TTree* treeDimuMC=((TTree*)fOutput->At(2));
  cout<<treeDimuMC->GetEntries()<<endl;
  TH1D *rapiditySpectrumMC1 = new TH1D("rapidity_distribution_2.5-4.0_MC","Rapidity distribution 2.5<#eta<4.0",3,2.5,4.0);
  // rapiditySpectrumMC1->Sumw2(kTRUE);
  treeDimuMC->Project("rapidity_distribution_2.5-4.0_MC","dimuon_rapidity","dimuon_rapidity>2.5 && dimuon_rapidity<4.0");
  rapiditySpectrumMC1->SetLineColor(kRed);
  rapiditySpectrumMC1->Draw("SAME E");

  canvHistos->cd(2);
  TH1D *ratio1 = (TH1D*)rapiditySpectrumMC1->Clone();
  ratio1->SetName("ratio1");
  ratio1->SetTitle("Ratio = #frac{Counts_{MC truth}-Counts_{Data}}{Counts_{Data}} 2.5<#eta<4.0");
  ratio1->GetYaxis()->SetTitle("Ratio");
  ratio1->Add(rapiditySpectrum1, -1.);
  ratio1->Divide(rapiditySpectrum1);
  ratio1->SetLineWidth(1);
  ratio1->SetLineColor(kBlack);
  ratio1->Draw("E");

  cout << "**********************" << endl;
  cout << "* Analysis completed *" << endl;
  cout << "**********************" << endl;

  return;
}
