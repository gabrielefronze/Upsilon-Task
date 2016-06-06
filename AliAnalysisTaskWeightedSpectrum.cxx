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
#include "TStyle.h"
#include "TLegend.h"
#include <vector>

using namespace std;
using namespace TMath;

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

  Double_t rapidityBinWidth = fRapidityAxis->GetBinWidth(1);
  Double_t rapidityLow = fRapidityAxis->GetXmin();

  if ( fMode==kFALSE ){

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
      static_cast<TF1*>(buffer)->SetName(Form("MCFittingFunction_%f-%f",rapidityLow+rapidityBinWidth*iInputResponseFunctions,rapidityLow+rapidityBinWidth*iInputResponseFunctions+rapidityBinWidth));
      if (!buffer) cout<<"problems"<<endl;
      fInputResponseFunctionsMC->Add(buffer);
      inputFile->GetObject(Form("dataFittingFunction_%d",iInputResponseFunctions),buffer);
      static_cast<TF1*>(buffer)->SetName(Form("dataFittingFunction_%f-%f",rapidityLow+rapidityBinWidth*iInputResponseFunctions,rapidityLow+rapidityBinWidth*iInputResponseFunctions+rapidityBinWidth));
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
  TH3D *histoMC = new TH3D("histo_MC","histo_MC",240,0.,12.,15,2.5,4.,320,0.,16.);
  TH3D *histoData = new TH3D("histo_data","histo_data",240,0.,12.,15,2.5,4.,320,0.,16.);

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
    Double_t rapidity1 = TMath::Abs(firstMuon->Eta());
    Double_t rapidity2 = TMath::Abs(secondMuon->Eta());
    Int_t rapidityBin1 = fRapidityAxis->FindBin(rapidity1) - 1;
    Int_t rapidityBin2 = fRapidityAxis->FindBin(rapidity2) - 1;

    // if ( rapidityBin1==10 ) rapidityBin1=9;
    // if ( rapidityBin2==10 ) rapidityBin2=9;

    //printf("%f, %f\n", rapidity1, rapidity2);

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

    if ( !fMode ){

      //cout<<"histo mode"<<endl;
      TH1D *weightingHisto1MC = static_cast<TH1D*>(weightingFunction1MC);
      TH1D *weightingHisto1Data = static_cast<TH1D*>(weightingFunction1Data);
      TH1D *weightingHisto2MC = static_cast<TH1D*>(weightingFunction2MC);
      TH1D *weightingHisto2Data = static_cast<TH1D*>(weightingFunction2Data);

      //cout<<weightingHisto1Data->GetName()<<endl;
      //cout<<weightingHisto2Data->GetName()<<endl;

      weightMC = weightingHisto1MC->GetBinContent(weightingHisto1MC->FindBin(pt1));
      weightData = weightingHisto1MC->GetBinContent(weightingHisto1Data->FindBin(pt1));
      weightMC *= weightingHisto2MC->GetBinContent(weightingHisto2MC->FindBin(pt2));
      weightData *= weightingHisto2MC->GetBinContent(weightingHisto2Data->FindBin(pt2));

    } else if ( fMode ){

      //cout<<"function mode"<<endl;

      TF1 *weightingFunctionMC1 = static_cast<TF1*>(weightingFunction1MC);
      TF1 *weightingFunctionData1 = static_cast<TF1*>(weightingFunction1Data);
      TF1 *weightingFunctionMC2 = static_cast<TF1*>(weightingFunction2MC);
      TF1 *weightingFunctionData2 = static_cast<TF1*>(weightingFunction2Data);

      //cout<<weightingFunctionData1->GetName()<<endl;
      //cout<<weightingFunctionData2->GetName()<<endl;

      weightMC = weightingFunctionMC1->Eval(pt1);
      weightData = weightingFunctionData1->Eval(pt1);

      weightMC *= weightingFunctionMC2->Eval(pt2);
      weightData *= weightingFunctionData2->Eval(pt2);
    }

    Double_t pt = dimuon.Pt();
    Double_t  rapidity = TMath::Abs(dimuon.Rapidity());
    Double_t mass = dimuon.M();

    histoMC->Fill(pt, rapidity, mass, weightMC);
    histoData->Fill(pt, rapidity, mass, weightData);

    //printf("%f, %f, %f, %f\n", pt, rapidity, mass, weightMC);
    //printf("%f, %f, %f, %f\n", pt, rapidity, mass, weightData);

  }

  delete muPlus;
  delete muMinus;

  PostData(1,fOutput);

  ////cout<<"########## ANALYSIS DONE! ##########"<<endl;

  return;
}

void AliAnalysisTaskWeightedSpectrum::Terminate(Option_t *) {

  gStyle->SetOptStat(0);
  gStyle->SetMarkerSize(0.8);

  fOutput=dynamic_cast<TList*>(GetOutputData(1));

  TH3D *histoMC=static_cast<TH3D*>(fOutput->At(0));
  TH3D *histoData=static_cast<TH3D*>(fOutput->At(1));
  histoMC->Sumw2(kTRUE);
  histoData->Sumw2(kTRUE);

  histoMC->GetXaxis()->SetTitle("P_{t} [Gev/c]");
  histoMC->GetYaxis()->SetTitle("Rapidity");
  histoMC->GetZaxis()->SetTitle("Mass [Gev/c^2]");
  histoData->GetXaxis()->SetTitle("P_{t} [Gev/c]");
  histoData->GetYaxis()->SetTitle("Rapidity");
  histoData->GetZaxis()->SetTitle("Mass [Gev/c^2]");
  //____________________________________________________________________________
  //histoMC->GetYaxis()->SetRangeUser(2.5,3.0);
  //histoData->GetYaxis()->SetRangeUser(2.5,3.0);

  TH1D *integratedSpectrumMC1 = histoMC->ProjectionY("integratedSpectrumMC1");
  integratedSpectrumMC1->SetName("MCSpectrum_2.5-3.0");
  integratedSpectrumMC1->SetTitle("Counts_{(MC,Data)} = Apt_{MC} #times #frac{Lpt_{(MC,Data)}}{Apt_{(MC,Data)}}  2.5<#eta<4.0");
  integratedSpectrumMC1->GetYaxis()->SetTitle("Counts");
  //integratedSpectrumMC1->SetTitle("2.5<#eta<3.0");
  integratedSpectrumMC1->SetLineColor(kRed);
  integratedSpectrumMC1->SetLineWidth(2);
  TH1D *integratedSpectrumData1 = histoData->ProjectionY("integratedSpectrumData1");
  integratedSpectrumData1->SetName("DataSpectrum_2.5-3.0");
  integratedSpectrumData1->SetTitle("DataSpectrum_2.5-3.0");
  integratedSpectrumMC1->Rebin(5);
  integratedSpectrumData1->Rebin(5);

  TH1D *ratio1 = (TH1D*)integratedSpectrumMC1->Clone();
  ratio1->SetName("ratio1");
  ratio1->SetTitle("Ratio = #frac{Counts_{MC}-Counts_{Data}}{Counts_{Data}} 2.5<#eta<4.0");
  ratio1->GetYaxis()->SetTitle("Ratio");
  ratio1->Add(integratedSpectrumData1, -1.);
  ratio1->Divide(integratedSpectrumData1);
  ratio1->SetLineWidth(1);
  ratio1->SetLineColor(kBlack);

  TLegend *leg1 = new TLegend(0.1,0.8,0.9,0.9);
  leg1->AddEntry(integratedSpectrumMC1,"Counts_{MC} (weight by simulated trigger response)","l");
  leg1->AddEntry(integratedSpectrumData1,"Counts_{Data} (weight by real trigger response)","l");

  TCanvas *canv1=new TCanvas("canv1","canv1");
  canv1->Divide(2);
  canv1->cd(1)->SetLogy();
  integratedSpectrumMC1->Draw("e");
  integratedSpectrumData1->Draw("same e");
  leg1->Draw();
  canv1->cd(2);
  ratio1->Draw("e");

  //____________________________________________________________________________
  //histoMC->GetYaxis()->SetRangeUser(3.0,3.5);
  //histoData->GetYaxis()->SetRangeUser(3.0,3.5);

  TH1D *integratedSpectrumMC2 = histoMC->ProjectionY("integratedSpectrumMC2");
  integratedSpectrumMC2->SetName("MCSpectrum_3.0-3.5");
  integratedSpectrumMC2->SetTitle("Counts_{(MC,Data)} = Apt_{MC} #times #frac{Lpt_{(MC,Data)}}{Apt_{(MC,Data)}}  2.5<#eta<4.0");
  integratedSpectrumMC2->GetYaxis()->SetTitle("Counts");
  //integratedSpectrumMC2->SetTitle("3.0<#eta<3.5");
  integratedSpectrumMC2->SetLineColor(kRed);
  integratedSpectrumMC2->SetLineWidth(2);
  TH1D *integratedSpectrumData2 = histoData->ProjectionY("integratedSpectrumData2");
  integratedSpectrumData2->SetName("DataSpectrum_3.0-3.5");
  integratedSpectrumData2->SetTitle("DataSpectrum_3.0-3.5");
  integratedSpectrumMC2->Rebin(5);
  integratedSpectrumData2->Rebin(5);

  TH1D *ratio2 = (TH1D*)integratedSpectrumMC2->Clone();
  ratio2->SetName("ratio2");
  ratio2->SetTitle("Ratio = #frac{Counts_{MC}-Counts_{Data}}{Counts_{Data}} 2.5<#eta<4.0");
  ratio2->GetYaxis()->SetTitle("Ratio");
  ratio2->Add(integratedSpectrumData2, -1.);
  ratio2->Divide(integratedSpectrumData2);
  ratio2->SetLineWidth(1);
  ratio2->SetLineColor(kBlack);

  TLegend *leg2 = new TLegend(0.1,0.8,0.9,0.9);
  leg2->AddEntry(integratedSpectrumMC2,"Counts_{MC} (weight by simulated trigger response)","l");
  leg2->AddEntry(integratedSpectrumData2,"Counts_{Data} (weight by real trigger response)","l");

  TCanvas *canv2=new TCanvas("canv2","canv2");
  canv2->Divide(2);
  canv2->cd(1)->SetLogy();
  integratedSpectrumMC2->Draw("e");
  integratedSpectrumData2->Draw("same e");
  leg2->Draw();
  canv2->cd(2);
  ratio2->Draw("e");

  //____________________________________________________________________________
  //histoMC->GetYaxis()->SetRangeUser(3.5,4.0);
  //histoData->GetYaxis()->SetRangeUser(3.5,4.0);

  TH1D *integratedSpectrumMC3 = histoMC->ProjectionY("integratedSpectrumMC3");
  integratedSpectrumMC3->SetName("MCSpectrum_3.5-4.0");
  integratedSpectrumMC3->SetTitle("Counts_{(MC,Data)} = Apt_{MC} #times #frac{Lpt_{(MC,Data)}}{Apt_{(MC,Data)}}  2.5<#eta<4.0");
  integratedSpectrumMC3->GetYaxis()->SetTitle("Counts");
  //integratedSpectrumMC3->SetTitle("3.5<#eta<4.0");
  integratedSpectrumMC3->SetLineColor(kRed);
  integratedSpectrumMC3->SetLineWidth(2);
  TH1D *integratedSpectrumData3 = histoData->ProjectionY("integratedSpectrumData3");
  integratedSpectrumData3->SetName("DataSpectrum_3.5-4.0");
  integratedSpectrumData3->SetTitle("DataSpectrum_3.5-4.0");
  integratedSpectrumMC3->Rebin(5);
  integratedSpectrumData3->Rebin(5);

  TH1D *ratio3 = (TH1D*)integratedSpectrumMC3->Clone();
  ratio3->SetName("ratio3");
  ratio3->SetTitle("Ratio = #frac{Counts_{MC}-Counts_{Data}}{Counts_{Data}} 2.5<#eta<4.0");
  ratio3->GetYaxis()->SetTitle("Ratio");
  ratio3->Add(integratedSpectrumData3, -1.);
  ratio3->Divide(integratedSpectrumData3);
  ratio3->SetLineWidth(1);
  ratio3->SetLineColor(kBlack);

  TLegend *leg3 = new TLegend(0.1,0.8,0.9,0.9);
  leg3->AddEntry(integratedSpectrumMC3,"Counts_{MC} (weight by simulated trigger response)","l");
  leg3->AddEntry(integratedSpectrumData3,"Counts_{Data} (weight by real trigger response)","l");

  TCanvas *canv3=new TCanvas("canv3","canv3");
  canv3->Divide(2);
  canv3->cd(1)->SetLogy();
  integratedSpectrumMC3->Draw("e");
  integratedSpectrumData3->Draw("same e");
  leg3->Draw();
  canv3->cd(2);
  ratio3->Draw("e");

  cout << "**********************" << endl;
  cout << "* Analysis completed *" << endl;
  cout << "**********************" << endl;

  return;
}
