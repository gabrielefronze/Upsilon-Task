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

void analysisMacroSingleCanv(){
	TFile *withPDCAFile=new TFile("full_analysis_outputs/withPDCA/AnalysisResults.root","READ");
	TFile *noPDCAFile=new TFile("full_analysis_outputs/noPDCA/AnalysisResults.root","READ");

	TFile *outputHistoFile=new TFile("outputHistos.root","RECREATE");


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

///////////////////////////////////////////////////////////////////

	TCanvas *canvHistos090=new TCanvas("canvHistos090","canvHistos090");
	canvHistos090->SetFrameBorderMode(0);
	TPad *Pad1=canvHistos090->cd();
	Pad1->Divide(1,2,0.01,0.0);
	Pad1->cd(1)->SetLogy();
	TH1D *histoInvariantMassWithPDCA090FullRap=new TH1D("Invariant mass 0%-90% 2.5<eta<4.0 PDCA","Invariant mass 0%-90% 2.5<eta<4.0 PDCA",250,2.5,15.);
	withPDCATree->Project("Invariant mass 0%-90% 2.5<eta<4.0 PDCA","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.");
//  histoInvariantMassWithPDCA090FullRap->ShowPeaks();
	histoInvariantMassWithPDCA090FullRap->Draw();
	outputHistoFile->cd();
	histoInvariantMassWithPDCA090FullRap->Write();

	TH1D *histoInvariantMassNoPDCA090FullRap=new TH1D("Invariant mass 0%-90% 2.5<eta<4.0","Invariant mass 0%-90% 2.5<eta<4.0",250,2.5,15.);
	noPDCATree->Project("Invariant mass 0%-90% 2.5<eta<4.0","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.");
//  histoInvariantMassNoPDCA090FullRap->ShowPeaks();
	histoInvariantMassNoPDCA090FullRap->SetLineColor(kRed);
	histoInvariantMassNoPDCA090FullRap->Draw("SAME");
	outputHistoFile->cd();
	histoInvariantMassNoPDCA090FullRap->Write();

	Pad1->cd(2);
	TH1D *histoInvariantMass090FullRapRatio=new TH1D("Ratio PDCA/NoPDCA 0%-90% 2.5<eta<4.0","Invariant mass 0%-90% 2.5<eta<4.0 PDCA",250,2.5,15.);
	histoInvariantMass090FullRapRatio->Add(histoInvariantMassWithPDCA090FullRap,histoInvariantMassNoPDCA090FullRap,1.,-1.*normalizationFactor);
	histoInvariantMass090FullRapRatio->Divide(histoInvariantMass090FullRapRatio,histoInvariantMassNoPDCA090FullRap);
	histoInvariantMass090FullRapRatio->SetLineColor(kBlack);
	histoInvariantMass090FullRapRatio->Draw();
	histoInvariantMass090FullRapRatio->Write();

///////////////////////////////////////////////////////////////////
	TCanvas *canvHistosAll=new TCanvas("canvHistosAll","canvHistosAll");
	canvHistosAll->SetFrameBorderMode(0);
	TPad *Pad2=canvHistosAll->cd();
	TH1D *histoInvariantMassWithPDCA=new TH1D("All entries PDCA","All entries PDCA",250,2.5,15.);
	withPDCATree->Project("All entries PDCA","dimuon_mass");
//  histoInvariantMassWithPDCA->ShowPeaks();
	histoInvariantMassWithPDCA->Draw();
	outputHistoFile->cd();
	histoInvariantMassWithPDCA->Write();

	TH1D *histoInvariantMassNoPDCA=new TH1D("All entries","All entries",250,2.5,15.);
	noPDCATree->Project("All entries","dimuon_mass");
//  histoInvariantMassNoPDCA->ShowPeaks();
	histoInvariantMassNoPDCA->SetLineColor(kRed);
	histoInvariantMassNoPDCA->Draw("SAME");  
	outputHistoFile->cd();
	histoInvariantMassNoPDCA->Write();  

///////////////////////////////////////////////////////////////////
	TCanvas *canvHistos020=new TCanvas("canvHistos020","canvHistos020");
	canvHistos020->SetFrameBorderMode(0);
	TPad *Pad3=canvHistos020->cd();
	Pad3->Divide(1,2,0.01,0.0);
	Pad3->cd(1)->SetLogy();
	TH1D *histoInvariantMassWithPDCA020FullRap=new TH1D("Invariant mass 0%-20% 2.5<eta<4.0 PDCA","Invariant mass 0%-20% 2.5<eta<4.0 PDCA",250,2.5,15.);
	withPDCATree->Project("Invariant mass 0%-20% 2.5<eta<4.0 PDCA","dimuon_mass","enevt_centrality>0 && enevt_centrality<20 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.");
//  histoInvariantMassWithPDCA020FullRap->ShowPeaks();
	histoInvariantMassWithPDCA020FullRap->Draw();
	outputHistoFile->cd();
	histoInvariantMassWithPDCA020FullRap->Write();

	TH1D *histoInvariantMassNoPDCA020FullRap=new TH1D("Invariant mass 0%-20% 2.5<eta<4.0","Invariant mass 0%-20% 2.5<eta<4.0",250,2.5,15.);
	noPDCATree->Project("Invariant mass 0%-20% 2.5<eta<4.0","dimuon_mass","enevt_centrality>0 && enevt_centrality<20 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.");
//  histoInvariantMassNoPDCA020FullRap->ShowPeaks();
	histoInvariantMassNoPDCA020FullRap->SetLineColor(kRed);
	histoInvariantMassNoPDCA020FullRap->Draw("SAME");  
	outputHistoFile->cd();
	histoInvariantMassNoPDCA020FullRap->Write();  

	Pad3->cd(2);
	TH1D *histoInvariantMass020FullRapRatio=new TH1D("Ratio PDCA/NoPDCA 0%-20% 2.5<eta<4.0","Invariant mass 0%-20% 2.5<eta<4.0 PDCA",250,2.5,15.);
	histoInvariantMass020FullRapRatio->Add(histoInvariantMassWithPDCA020FullRap,histoInvariantMassNoPDCA020FullRap,1.,-1.*normalizationFactor);
	histoInvariantMass020FullRapRatio->Divide(histoInvariantMass020FullRapRatio,histoInvariantMassNoPDCA020FullRap);
	histoInvariantMass020FullRapRatio->SetLineColor(kBlack);
	histoInvariantMass020FullRapRatio->Draw();
	outputHistoFile->cd();
	histoInvariantMass020FullRapRatio->Write();

///////////////////////////////////////////////////////////////////
	TCanvas *canvHistos2090=new TCanvas("canvHistos2090","canvHistos2090");
	canvHistos2090->SetFrameBorderMode(0);
	TPad *Pad4=canvHistos2090->cd();
	Pad4->Divide(1,2,0.01,0.0);
	Pad4->cd(1)->SetLogy();
	TH1D *histoInvariantMassWithPDCA2090FullRap=new TH1D("20%-90% 2.5<eta<4.0 PDCA","20%-90% 2.5<eta<4.0 PDCA",250,2.5,15.);
	withPDCATree->Project("20%-90% 2.5<eta<4.0 PDCA","dimuon_mass","enevt_centrality>20 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.");
//  histoInvariantMassWithPDCA2090FullRap->ShowPeaks();
	histoInvariantMassWithPDCA2090FullRap->Draw();
	outputHistoFile->cd();
	histoInvariantMassWithPDCA2090FullRap->Write();

	TH1D *histoInvariantMassNoPDCA2090FullRap=new TH1D("20%-90% 2.5<eta<4.0","20%-90% 2.5<eta<4.0",250,2.5,15.);
	noPDCATree->Project("20%-90% 2.5<eta<4.0","dimuon_mass","enevt_centrality>20 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.");
//  histoInvariantMassNoPDCA2090FullRap->ShowPeaks();
	histoInvariantMassNoPDCA2090FullRap->SetLineColor(kRed);
	histoInvariantMassNoPDCA2090FullRap->Draw("SAME");  
	outputHistoFile->cd();
	histoInvariantMassNoPDCA2090FullRap->Write();  

	Pad4->cd(2);
	TH1D *histoInvariantMass2090FullRapRatio=new TH1D("20%-90% 2.5<eta<4.0","20%-90% 2.5<eta<4.0 PDCA",250,2.5,15.);
	histoInvariantMass2090FullRapRatio->Add(histoInvariantMassWithPDCA2090FullRap,histoInvariantMassNoPDCA2090FullRap,1.,-1.*normalizationFactor);
	histoInvariantMass2090FullRapRatio->Divide(histoInvariantMass2090FullRapRatio,histoInvariantMassNoPDCA2090FullRap);
	histoInvariantMass2090FullRapRatio->SetLineColor(kBlack);
	histoInvariantMass2090FullRapRatio->Draw();
	outputHistoFile->cd();
	histoInvariantMass2090FullRapRatio->Write();

///////////////////////////////////////////////////////////////////
	TCanvas *canvHistos090r1=new TCanvas("canvHistos090r1","canvHistos090r1");
	canvHistos090r1->SetFrameBorderMode(0);
	TPad *Pad5=canvHistos090r1->cd();
	Pad5->Divide(1,2,0.01,0.0);
	Pad5->cd(1)->SetLogy();
	TH1D *histoInvariantMassWithPDCA090SecondRap=new TH1D("Invariant mass 0%-90% 2.5<eta<3.2 PDCA","Invariant mass 0%-90% 2.5<eta<3.2 PDCA",250,2.5,15.);
	withPDCATree->Project("Invariant mass 0%-90% 2.5<eta<3.2 PDCA","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-3.2 && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.");
//  histoInvariantMassWithPDCA090SecondRap->ShowPeaks();
	histoInvariantMassWithPDCA090SecondRap->Draw();
	outputHistoFile->cd();
	histoInvariantMassWithPDCA090SecondRap->Write();

	TH1D *histoInvariantMassNoPDCA090SecondRap=new TH1D("Invariant mass 0%-90% 2.5<eta<3.2","Invariant mass 0%-90% 2.5<eta<3.2",250,2.5,15.);
	noPDCATree->Project("Invariant mass 0%-90% 2.5<eta<3.2","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-3.2 && dimuon_rapidity<-2.5 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.");
//  histoInvariantMassNoPDCA090SecondRap->ShowPeaks();
	histoInvariantMassNoPDCA090SecondRap->SetLineColor(kRed);
	histoInvariantMassNoPDCA090SecondRap->Draw("SAME");
	outputHistoFile->cd();
	histoInvariantMassNoPDCA090SecondRap->Write();

	Pad5->cd(2);
	TH1D *histoInvariantMass090SecondRapRatio=new TH1D("Ratio PDCA/NoPDCA 0%-90% 2.5<eta<3.2","Invariant mass 0%-90% 2.5<eta<3.2 PDCA",250,2.5,15.);
	histoInvariantMass090SecondRapRatio->Add(histoInvariantMassWithPDCA090SecondRap,histoInvariantMassNoPDCA090SecondRap,1.,-1.*normalizationFactor);
	histoInvariantMass090SecondRapRatio->Divide(histoInvariantMass090SecondRapRatio,histoInvariantMassNoPDCA090SecondRap);
	histoInvariantMass090SecondRapRatio->SetLineColor(kBlack);
	histoInvariantMass090SecondRapRatio->Draw();
	outputHistoFile->cd();
	histoInvariantMass090SecondRapRatio->Write();

///////////////////////////////////////////////////////////////////
	TCanvas *canvHistos090r2=new TCanvas("canvHistos090r2","canvHistos090r2");
	canvHistos090r2->SetFrameBorderMode(0);
	TPad *Pad6=canvHistos090r2->cd();
	Pad6->Divide(1,2,0.01,0.0);
	Pad6->cd(1)->SetLogy();
	TH1D *histoInvariantMassWithPDCA090FirstRap=new TH1D("Invariant mass 0%-90% 3.2<eta<4.0 PDCA","Invariant mass 0%-90% 3.2<eta<4.0 PDCA",250,2.5,15.);
	withPDCATree->Project("Invariant mass 0%-90% 3.2<eta<4.0 PDCA","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-3.2 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.");
//  histoInvariantMassWithPDCA090FirstRap->ShowPeaks();
	histoInvariantMassWithPDCA090FirstRap->Draw();
	outputHistoFile->cd();
	histoInvariantMassWithPDCA090FirstRap->Write();

	TH1D *histoInvariantMassNoPDCA090FirstRap=new TH1D("Invariant mass 0%-90% 3.2<eta<4.0","Invariant mass 0%-90% 3.2<eta<4.0",250,2.5,15.);
	noPDCATree->Project("Invariant mass 0%-90% 3.2<eta<4.0","dimuon_mass","enevt_centrality>0 && enevt_centrality<90 && dimuon_rapidity>-4. && dimuon_rapidity<-3.2 && highest_muon_transverse_momentum>0. && lowest_muon_transverse_momentum>0.");
//  histoInvariantMassNoPDCA090FirstRap->ShowPeaks();
	histoInvariantMassNoPDCA090FirstRap->SetLineColor(kRed);
	histoInvariantMassNoPDCA090FirstRap->Draw("SAME");
	outputHistoFile->cd();
	histoInvariantMassNoPDCA090FirstRap->Write();

	Pad6->cd(2);
	TH1D *histoInvariantMass090FirstRapRatio=new TH1D("Ratio PDCA/NoPDCA 0%-90% 3.2<eta<4.0","Invariant mass 0%-90% 3.2<eta<4.0 PDCA",250,2.5,15.);
	histoInvariantMass090FirstRapRatio->Add(histoInvariantMassWithPDCA090FirstRap,histoInvariantMassNoPDCA090FirstRap,1.,-1.*normalizationFactor);
	histoInvariantMass090FirstRapRatio->Divide(histoInvariantMass090FirstRapRatio,histoInvariantMassNoPDCA090FirstRap);
	histoInvariantMass090FirstRapRatio->SetLineColor(kBlack);
	histoInvariantMass090FirstRapRatio->Draw();
	outputHistoFile->cd();
	histoInvariantMass090FirstRapRatio->Write();
}