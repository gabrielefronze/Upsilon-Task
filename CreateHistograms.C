#include "TH1.h"
#include "TMath.h"
#include "TTree.h"
#include "TList.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include <fstream>
#include <iostream>

using namespace std;

const Int_t MAX_CHARS_PER_LINE=512;

Double_t *fTreeData;
enum{
    kMomentum,
    kTransverse,
    kRapidity,
    kCentrality,
    kLowMomentumMuonpt,
    kHighMomentumMuonpt,
    kMass,
    kSparseDimension
};


TCanvas *canvHistos=new TCanvas("canvHistos","canvHistos");

TH1D *histoInvariantMass090FullRap=new TH1D("0%-90% 2.5<eta<4.0","0%-90% 2.5<eta<4.0",340,2.,12.);
TH1D *histoInvariantMass020FullRap=new TH1D("0%-20% 2.5<eta<4.0","0%-20% 2.5<eta<4.0",340,2.,12.);
TH1D *histoInvariantMass2090FullRap=new TH1D("20%-90% 2.5<eta<4.0","20%-90% 2.5<eta<4.0",340,2.,12.);
TH1D *histoInvariantMass090SecondRap=new TH1D("0%-90% 2.5<eta<3.2","0%-90% 2.5<eta<3.2",340,2.,12.);
TH1D *histoInvariantMass090FirstRap=new TH1D("0%-90% 3.2<eta<4.0","0%-90% 3.2<eta<4.0",340,2.,12.);

void LoopOverEntries(TFile *inputfile){
	TList *list;
	TTree *tree;
	fTreeData=new Double_t[kSparseDimension];

	inputfile->GetObject("Upsilonfirst/UpsilonOutfirst",list);
	if(!list) return;

	tree=(TTree*)list->FindObject("MuonData");
	Int_t entries=tree->GetEntries();
	if(!tree || entries==0) return;

	Bool_t control=kTRUE;
	TBranch *branchBuffer;
	branchBuffer=tree->GetBranch("dimuon_momentum");
	if(!branchBuffer) control=kFALSE;
	branchBuffer->SetAddress(&fTreeData[kMomentum]);

	branchBuffer=tree->GetBranch("dimuon_transverse_momentum");
	if(!branchBuffer) control=kFALSE;
	branchBuffer->SetAddress(&fTreeData[kTransverse]);

	branchBuffer=tree->GetBranch("dimuon_rapidity");
	if(!branchBuffer) control=kFALSE;
	branchBuffer->SetAddress(&fTreeData[kRapidity]);

	branchBuffer=tree->GetBranch("dimuon_mass");
	if(!branchBuffer) control=kFALSE;
	branchBuffer->SetAddress(&fTreeData[kMass]);

	branchBuffer=tree->GetBranch("highest_muon_transverse_momentum");
	if(!branchBuffer) control=kFALSE;
	branchBuffer->SetAddress(&fTreeData[kHighMomentumMuonpt]);

	branchBuffer=tree->GetBranch("lowest_muon_transverse_momentum");
	if(!branchBuffer) control=kFALSE;
	branchBuffer->SetAddress(&fTreeData[kLowMomentumMuonpt]);

	branchBuffer=tree->GetBranch("enevt_centrality");
	if(!branchBuffer) control=kFALSE;
	branchBuffer->SetAddress(&fTreeData[kCentrality]);

	if(!control) return;

	for(Int_t i=0;i<entries;i++){
		tree->GetEvent(i);

		histoInvariantMass090FullRap->Fill(fTreeData[kMass]);
		
		//if(i%10000==0) cout<<i<<"/"<<entries<<endl;

		if(fTreeData[kCentrality]>0. && fTreeData[kCentrality]<90.){
			if(fTreeData[kRapidity]>-4. && fTreeData[kRapidity]<-2.5){
				//histoInvariantMass090FullRap->Fill(fTreeData[kMass]);
			}
			if(fTreeData[kRapidity]>-4. && fTreeData[kRapidity]<-3.2){
				histoInvariantMass090FirstRap->Fill(fTreeData[kMass]);
			}
			if(fTreeData[kRapidity]>-3.2 && fTreeData[kRapidity]<-2.5){
				histoInvariantMass090SecondRap->Fill(fTreeData[kMass]);
			}
		}

		if(fTreeData[kCentrality]>20. && fTreeData[kCentrality]<90.){
			if(fTreeData[kRapidity]>-4. && fTreeData[kRapidity]<-2.5){
				histoInvariantMass2090FullRap->Fill(fTreeData[kMass]);
			}
		}

		if(fTreeData[kCentrality]>0. && fTreeData[kCentrality]<20.){
			if(fTreeData[kRapidity]>-4. && fTreeData[kRapidity]<-2.5){
				histoInvariantMass020FullRap->Fill(fTreeData[kMass]);
			}
		}
	}

	delete[] fTreeData;

	canvHistos->cd(1);
	histoInvariantMass090FullRap->Draw();

	canvHistos->cd(3);
	histoInvariantMass020FullRap->Draw();

	canvHistos->cd(4);
	histoInvariantMass2090FullRap->Draw();

	canvHistos->cd(5);
	histoInvariantMass090SecondRap->Draw();

	canvHistos->cd(6);
	histoInvariantMass090FirstRap->Draw();

	canvHistos->Update();
	canvHistos->Modified();

}

void CreateHistograms(TString fileList){
	ifstream fin;

	canvHistos->Divide(2,3);

    TFile *inputfile;
    fin.open(fileList.Data());
    while(!fin.eof()){
    	string filename;
        getline(fin,filename);
        //cout<<"filename "<<filename<<endl;
        //if(fin.eof())break;
        inputfile=new TFile(filename.c_str());
        if(!inputfile) continue;
        cout<<"filename "<<filename<<endl;
        LoopOverEntries(inputfile);
        inputfile=0x0;
    }
    Printf("\t\t Run adding ended.\n");
    fin.close();
}