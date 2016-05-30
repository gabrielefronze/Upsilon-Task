#include "TH1.h"
#include "TMath.h"
#include "TTree.h"
#include "TList.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include <fstream>
#include <iostream>

void analysisMacro2(){
	TFile *withPDCAFile=new TFile("full_analysis_outputs/withPDCA175/AnalysisResults.root","READ");
	TFile *noPDCAFile=new TFile("full_analysis_outputs/noPDCA/AnalysisResults.root","READ");

	TFile *outputHistoFile=new TFile("outputHistosLardeaux.root","RECREATE");


	TList *withPDCAProcessedNumberList;
	withPDCAFile->GetObject("MultSelection/cListMultSelection",withPDCAProcessedNumberList);
	TList *withPDCAList;
	withPDCAFile->GetObject("Upsilonfirst/UpsilonOutfirst",withPDCAList);

	TList *noPDCAProcessedNumberList;
	noPDCAFile->GetObject("MultSelection/cListMultSelection",noPDCAProcessedNumberList);
	TList *noPDCAList;
	noPDCAFile->GetObject("Upsilonfirst/UpsilonOutfirst",noPDCAList);


	TH1D *withPDCAProcessedHisto;
	withPDCAProcessedHisto=(TH1D*)withPDCAProcessedNumberList->FindObject("fHistEventCounter");
	TTree *withPDCATree;
	withPDCATree=(TTree*)withPDCAList->FindObject("MuonData");

	TH1D *noPDCAProcessedHisto;
	noPDCAProcessedHisto=(TH1D*)noPDCAProcessedNumberList->FindObject("fHistEventCounter");
	TTree *noPDCATree;
	noPDCATree=(TTree*)noPDCAList->FindObject("MuonData");

	Double_t processedWithPDCA=withPDCAProcessedHisto->GetBinContent(1);
	Double_t processedNoPDCA=noPDCAProcessedHisto->GetBinContent(1);
	Double_t normalizationFactor=processedWithPDCA/processedNoPDCA;

	Int_t centralityBins[5]={0,10,30,50,90};

	for(Int_t iCentrality=0; iCentrality<4; iCentrality++){
		///////////////////////////////////////////////////////////////////
		TH1D *histoPDCA=new TH1D(Form("%d-%d 2.5<eta<4.0 PDCA",centralityBins[iCentrality],centralityBins[iCentrality+1]),Form("%d-%d 2.5<eta<4.0 PDCA",centralityBins[iCentrality],centralityBins[iCentrality+1]),250,2.5,15.);
		withPDCATree->Project(Form("%d-%d 2.5<eta<4.0 PDCA",centralityBins[iCentrality],centralityBins[iCentrality+1]),"dimuon_mass",Form("enevt_centrality>%d && enevt_centrality<%d && dimuon_rapidity>-4. && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.",centralityBins[iCentrality],centralityBins[iCentrality+1]));
		outputHistoFile->cd();
		histoPDCA->Write();

		TH1D *histoNoPDCA=new TH1D(Form("%d-%d 2.5<eta<4.0",centralityBins[iCentrality],centralityBins[iCentrality+1]),Form("%d-%d 2.5<eta<4.0",centralityBins[iCentrality],centralityBins[iCentrality+1]),250,2.5,15.);
		noPDCATree->Project(Form("%d-%d 2.5<eta<4.0",centralityBins[iCentrality],centralityBins[iCentrality+1]),"dimuon_mass",Form("enevt_centrality>%d && enevt_centrality<%d && dimuon_rapidity>-4. && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.",centralityBins[iCentrality],centralityBins[iCentrality+1]));
		histoNoPDCA->SetLineColor(kRed);
		outputHistoFile->cd();
		histoNoPDCA->Write();

		TH1D *histoRatio=new TH1D(Form("%d-%d 2.5<eta<4.0 ratio",centralityBins[iCentrality],centralityBins[iCentrality+1]),Form("%d-%d 2.5<eta<4.0 ratio",centralityBins[iCentrality],centralityBins[iCentrality+1]),250,2.5,15.);
		histoRatio->Add(histoPDCA,histoNoPDCA,1.,-1.*normalizationFactor);
		histoRatio->Divide(histoRatio,histoNoPDCA);
		histoRatio->SetLineColor(kBlack);
		histoRatio->Write();

		TCanvas *canv=new TCanvas(Form("canvas %d-%d 2.5<eta<4.0",centralityBins[iCentrality],centralityBins[iCentrality+1]),Form("canvas %d-%d 2.5<eta<4.0",centralityBins[iCentrality],centralityBins[iCentrality+1]));
		canv->Divide(1,2);
		canv->cd(1)->SetLogy();
		histoPDCA->DrawClone();
		histoNoPDCA->DrawClone("SAME");
		canv->cd(2);
		histoRatio->DrawClone();
		outputHistoFile->cd();
		canv->Write();

		delete histoPDCA;
		delete histoNoPDCA;
		delete histoRatio;
		delete canv;
	}

	for(Int_t iCentrality=1; iCentrality<4; iCentrality++){
		///////////////////////////////////////////////////////////////////
		TH1D *histoPDCA=new TH1D(Form("%d-%d 2.5<eta<4.0 PDCA",centralityBins[0],centralityBins[iCentrality+1]),Form("%d-%d 2.5<eta<4.0 PDCA",centralityBins[0],centralityBins[iCentrality+1]),250,2.5,15.);
		withPDCATree->Project(Form("%d-%d 2.5<eta<4.0 PDCA",centralityBins[0],centralityBins[iCentrality+1]),"dimuon_mass",Form("enevt_centrality>%d && enevt_centrality<%d && dimuon_rapidity>-4. && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.",centralityBins[0],centralityBins[iCentrality+1]));
		outputHistoFile->cd();
		histoPDCA->Write();

		TH1D *histoNoPDCA=new TH1D(Form("%d-%d 2.5<eta<4.0",centralityBins[0],centralityBins[iCentrality+1]),Form("%d-%d 2.5<eta<4.0",centralityBins[0],centralityBins[iCentrality+1]),250,2.5,15.);
		noPDCATree->Project(Form("%d-%d 2.5<eta<4.0",centralityBins[0],centralityBins[iCentrality+1]),"dimuon_mass",Form("enevt_centrality>%d && enevt_centrality<%d && dimuon_rapidity>-4. && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.",centralityBins[0],centralityBins[iCentrality+1]));
		histoNoPDCA->SetLineColor(kRed);
		outputHistoFile->cd();
		histoNoPDCA->Write();

		TH1D *histoRatio=new TH1D(Form("%d-%d 2.5<eta<4.0 ratio",centralityBins[0],centralityBins[iCentrality+1]),Form("%d-%d 2.5<eta<4.0 ratio",centralityBins[0],centralityBins[iCentrality+1]),250,2.5,15.);
		histoRatio->Add(histoPDCA,histoNoPDCA,1.,-1.*normalizationFactor);
		histoRatio->Divide(histoRatio,histoNoPDCA);
		histoRatio->SetLineColor(kBlack);
		histoRatio->Write();

		TCanvas *canv=new TCanvas(Form("canvas %d-%d 2.5<eta<4.0",centralityBins[0],centralityBins[iCentrality+1]),Form("canvas %d-%d 2.5<eta<4.0",centralityBins[0],centralityBins[iCentrality+1]));
		canv->Divide(1,2);
		canv->cd(1)->SetLogy();
		histoPDCA->DrawClone();
		histoNoPDCA->DrawClone("SAME");
		canv->cd(2);
		histoRatio->DrawClone();
		outputHistoFile->cd();
		canv->Write();

		delete histoPDCA;
		delete histoNoPDCA;
		delete histoRatio;
		delete canv;
	}
}