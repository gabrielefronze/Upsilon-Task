#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "iostream"
#include "TStyle.h"

using namespace std;

void AxEffSyst(TString modOCDBOutputPath, TString stdOCDBOutputPath){

  gStyle->SetOptStat(0);

  TFile *axEffModOCDBFile = new TFile(modOCDBOutputPath.Data(),"READONLY");
  TFile *axEffStdOCDBFile = new TFile(stdOCDBOutputPath.Data(),"READONLY");

  TH1D *axEffModOCDBHisto = 0x0;
  axEffModOCDBFile->GetObject("ratio1", axEffModOCDBHisto);
  axEffModOCDBHisto->SetLineColor(kBlue);
  axEffModOCDBHisto->SetLineWidth(2);
  axEffModOCDBHisto->SetTitle("A#times#epsilon");
  axEffModOCDBHisto->GetXaxis()->SetTitle("Rapidity");
  axEffModOCDBHisto->GetYaxis()->SetTitle("A#times#epsilon");

  TH1D *axEffStdOCDBHisto = 0x0;
  axEffStdOCDBFile->GetObject("ratio1", axEffStdOCDBHisto);
  axEffStdOCDBHisto->SetLineColor(kRed);
  axEffStdOCDBHisto->SetTitle("A#times#epsilon Std OCDB");
  axEffStdOCDBHisto->GetXaxis()->SetTitle("Rapidity");
  axEffStdOCDBHisto->GetYaxis()->SetTitle("A#times#epsilon");

  axEffModOCDBHisto->Sumw2(kTRUE);
  axEffStdOCDBHisto->Sumw2(kTRUE);

  TH1D *ratio = (TH1D*)axEffStdOCDBHisto->Clone();
  ratio->Sumw2(kTRUE);
  ratio->SetLineWidth(1);
  ratio->Add(axEffModOCDBHisto, -1.);
  ratio->Divide(axEffStdOCDBHisto);
  ratio->SetLineColor(kBlack);
  ratio->SetTitle("#frac{A#times#epsilon_{Std OCDB}-A#times#epsilon_{Mod OCDB}}{A#times#epsilon_{Std OCDB}}");
  ratio->GetXaxis()->SetTitle("Rapidity");
  ratio->GetYaxis()->SetTitle("Ratio (%)");

  for (Int_t iBinsRatio = 0; iBinsRatio <= ratio->GetNbinsX()+1; iBinsRatio++) {
    Double_t binContent = ratio->GetBinContent(iBinsRatio);
    cout<<binContent<<"->";
    if ( binContent<0. ) ratio->SetBinContent(iBinsRatio, -binContent*100.);
    else ratio->SetBinContent(iBinsRatio, binContent*100.);
    binContent = ratio->GetBinContent(iBinsRatio);
    cout<<binContent<<endl;
  }

  TCanvas *canv = new TCanvas("canv");
  canv->Divide(2,1);

  canv->cd(1);
  axEffModOCDBHisto->SetDirectory(0);
  axEffModOCDBHisto->Draw("E");
  axEffStdOCDBHisto->SetDirectory(0);
  axEffStdOCDBHisto->Draw("SAME E");

  canv->cd(2);
  ratio->GetYaxis()->SetRangeUser(-0.5,ratio->GetMaximum()*2);
  ratio->SetDirectory(0);
  ratio->Draw("E");

  axEffModOCDBFile->Close();
  axEffStdOCDBFile->Close();
}
