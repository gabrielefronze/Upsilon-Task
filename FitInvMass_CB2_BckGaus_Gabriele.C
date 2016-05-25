#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TASImage.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMatrixDSym.h>
#include <TFitResultPtr.h>
#endif

Double_t FuncJpsi(Double_t *, Double_t *);
Double_t FuncJpsi2(Double_t *, Double_t *);
Double_t FuncPsiP(Double_t *, Double_t *);
Double_t FuncPsiP2(Double_t *, Double_t *);
Double_t FuncBck(Double_t *, Double_t *);
Double_t FuncTot(Double_t *, Double_t *);

Double_t NormPsiP=0;

void FitInvMass_CB2_BckGaus_Gabriele(TString inputFile, TString histoName, char *period = "LHC15o", char *tailsfix="yes"){

//-----------------------------------------------
// style
//-----------------------------------------------
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
gStyle->SetCanvasColor(10);
gStyle->SetFrameFillColor(10);

//--------------------------------------------
// Open file
//--------------------------------------------
TFile *f = new TFile(inputFile.Data());
printf("Opening %s\n",inputFile.Data());

TH1D *histo;
histo = new TH1D("histo","histo",250,2.5,15.);
Int_t ptbin;
histo->Add((TH1D*) f->Get(histoName.Data()));  
printf("Using Histo= %s\n",histoName.Data());

histo->Draw();

//-----------------------------------------------------------------------
// choose the JPsi tails
//-----------------------------------------------------------------------

Double_t cbalpha1; 
Double_t cbn1;
Double_t cbalpha2;
Double_t cbn2;

cbalpha1=1.03; 
cbn1=4.36;
cbalpha2=2.22;
cbn2=3.17;  

printf("\nDouble CB tails:\n");   
printf("\nAlpha1=%f N1=%f Alpha2=%f N2=%f\n",cbalpha1,cbn1,cbalpha2,cbn2);

//-----------------------------------------------------------------------
// prepare plot
//-----------------------------------------------------------------------
//histo->Rebin(2);

histo->GetXaxis()->SetTitle("M_{#mu#mu} (GeV/c^{2})");
histo->GetYaxis()->SetTitle("dN/dM_{#mu#mu}");
Double_t binWidth= histo->GetBinWidth(1);

TCanvas *c = new TCanvas("c","c",20,20,600,600);
c->cd(1);
// set to 0 the entries of the JPsi region (2.9-3.3)
TH1D *histoForBck = (TH1D*) histo->Clone("histoForBck");
Double_t xxmin=2.9/histo->GetBinWidth(1);
Double_t xxmax=3.3/histo->GetBinWidth(1);

for(int i=xxmin;i<xxmax;i++){
  histoForBck->SetBinContent(i+1,0);
  histoForBck->SetBinError(i+1,0);
}
//---------------------------------------------------------------------
// 1st step: Bck fit
//---------------------------------------------------------------------
Double_t par[12]={90000.,1.2,0.6,0.2,5000.,3.096,7.0e-02,1.00011e+00,3.70078e+00,1.68359e+00,3.63002e+00,0.01}; 
TF1 *funcbck;
funcbck = new TF1("funcbck",FuncBck,2.2,5.5,4); 
funcbck->SetParameter(0,par[0]);
funcbck->SetParameter(1,par[1]);
funcbck->SetParameter(2,par[2]);
funcbck->SetParameter(3,par[3]);
funcbck->SetLineColor(kOrange);
funcbck->SetLineWidth(2);
printf("\n--- Step 3: Fitting the bck1+bck2 ---\n");
histoForBck->Fit(funcbck,"R");

//---------------------------------------------------------------------
// 2nd step: JPsi fit
//---------------------------------------------------------------------

TF1 *funcpsi = new TF1("funcpsi",FuncJpsi,2.9,3.3,12);
funcpsi->SetParameters(par);
funcpsi->FixParameter(0,funcbck->GetParameter(0));
funcpsi->FixParameter(1,funcbck->GetParameter(1));
funcpsi->FixParameter(2,funcbck->GetParameter(2));
funcpsi->FixParameter(3,funcbck->GetParameter(3));
funcpsi->SetParameter(5,par[5]);
funcpsi->FixParameter(6,par[6]);
funcpsi->FixParameter(7,cbalpha1);
funcpsi->FixParameter(8,cbn1);
funcpsi->FixParameter(9,cbalpha2);
funcpsi->FixParameter(10,cbn2);
funcpsi->FixParameter(11,par[11]);
funcpsi->SetParLimits(4,0.,100.);  //fixing sigma defined positive!

funcpsi->SetLineColor(kGreen);
funcpsi->SetLineWidth(2);
printf("\n--- Step 1: Fitting the JPsi ---\n");
histo->Fit(funcpsi,"R"); 

//---------------------------------------------------------------------
// 3rd step: total fit
//---------------------------------------------------------------------
Double_t FitMin=2.2, FitMax=4.5;
TF1 *functot = new TF1("functot",FuncTot,FitMin,FitMax,12); 
functot->SetParameter(0,funcbck->GetParameter(0));   
functot->SetParameter(1,funcbck->GetParameter(1));   
functot->SetParameter(2,funcbck->GetParameter(2));   
functot->SetParameter(3,funcbck->GetParameter(3));   
functot->SetParameter(4,funcpsi->GetParameter(4));   
functot->SetParameter(5,3.096); 
functot->SetParameter(6,0.07);

if((strcmp(tailsfix,"yes")==0)){
  functot->FixParameter(7,cbalpha1); 
  functot->FixParameter(8,cbn1);     
  functot->FixParameter(9,cbalpha2); 
  functot->FixParameter(10,cbn2);     
}  else if((strcmp(tailsfix,"no")==0)){
  functot->SetParameter(7,cbalpha1); 
  functot->SetParameter(8,cbn1);     
  functot->SetParameter(9,cbalpha2); 
  functot->SetParameter(10,cbn2);     
  functot->SetParLimits(7,0.,cbalpha1*10); 
  functot->SetParLimits(8,0.,cbn1*10);     
  functot->SetParLimits(9,0.,cbalpha2*10); 
  functot->SetParLimits(10,0.,cbn2*10);     
}

functot->SetParameter(11,0.01);   
functot->SetParLimits(11,0.,1.);   

functot->SetLineColor(kBlue);
functot->SetLineWidth(2);

printf("\n--- Step 4: Fitting total spectrum ---\n");
TFitResultPtr r = histo->Fit(functot,"RLS");
TMatrixDSym cov = r->GetCovarianceMatrix();
cov.Print();

TMatrixDSym cov2 = r->GetCorrelationMatrix();
cov2.Print();

Double_t *fullmat;
fullmat = cov.GetMatrixArray();
Double_t psimat[49];
for(Int_t i=0;i<7;i++){
  for(Int_t j=0;j<7;j++){
    psimat[7*i+j]=fullmat[52+j+12*i];
  }
}

histo->GetXaxis()->SetRangeUser(2.,5.);
histo->SetMinimum(1.);
gPad->SetLogy(1);
histo->Draw("e");
histo->SetMarkerStyle(20);
histo->SetMarkerColor(2);
histo->SetMarkerSize(0.7);
histo->SetMinimum(1);
Double_t chi2 = functot->GetChisquare();
Double_t ndf = functot->GetNDF();
printf("Chi2/ndf= %f\n",chi2/ndf);

//------------------------------------------------------
// plot fitting functions and compute uncertainties
//------------------------------------------------------
TF1 *psifix = new TF1("psifix",FuncJpsi,0.,5.,12); 
psifix->SetParameter(0,functot->GetParameter(0));
psifix->SetParameter(1,functot->GetParameter(1));
psifix->SetParameter(2,functot->GetParameter(2));
psifix->SetParameter(3,functot->GetParameter(3)); 
psifix->SetParameter(4,functot->GetParameter(4));  
psifix->SetParameter(5,functot->GetParameter(5));  
psifix->SetParameter(6,functot->GetParameter(6));  
psifix->SetParameter(7,functot->GetParameter(7));  
psifix->SetParameter(8,functot->GetParameter(8));  
psifix->SetParameter(9,functot->GetParameter(9));  
psifix->SetParameter(10,functot->GetParameter(10));  
psifix->SetParameter(11,functot->GetParameter(11));  
psifix->SetLineColor(kGreen);
psifix->Draw("same");
Double_t NPsi=psifix->Integral(0.,5.)/binWidth;

TF1 *psifix2 = new TF1("psifix2",FuncJpsi2,0.,5.,7); 
psifix2->SetParameter(0,functot->GetParameter(4));
psifix2->SetParameter(1,functot->GetParameter(5));
psifix2->SetParameter(2,functot->GetParameter(6));
psifix2->SetParameter(3,functot->GetParameter(7)); 
psifix2->SetParameter(4,functot->GetParameter(8));  
psifix2->SetParameter(5,functot->GetParameter(9));  
psifix2->SetParameter(6,functot->GetParameter(10));  
Double_t psipar[7];
for(int i=0;i<7;i++) {
  psipar[i]=functot->GetParameter(4+i);
} 
Double_t ErrPsiCorrParam = psifix2->IntegralError(0.,5.,psipar,psimat)/binWidth;
Double_t ErrNPsi=ErrPsiCorrParam;
printf("\nPsi=%f \n",psifix2->Integral(0.,5.)/binWidth);
printf("ErrPsi=%4.3f (computed using corr. params) %f, %f\n",ErrPsiCorrParam, psifix2->IntegralError(0.,5.,psipar,psimat),binWidth);

char text3[100];
sprintf(text3,"N_{J/#psi}= %3.0f #pm %2.0f",NPsi,ErrPsiCorrParam);
TLatex *l1 = new TLatex(3.5,psifix2->Integral(0.,5.)*0.25/binWidth,text3);
l1->SetTextColor(kBlue);
l1->SetTextSize(0.038);
l1->SetTextFont(42);
l1->Draw();

char text4[100];
sprintf(text4,"m_{J/#psi}= %5.3f #pm %5.3f GeV/c^{2}",functot->GetParameter(5),functot->GetParError(5));
TLatex *l2 = new TLatex(3.5,psifix2->Integral(0.,5.)/binWidth*0.1,text4);
l2->SetTextColor(kBlue);
l2->SetTextSize(0.035);
l2->SetTextFont(42);
l2->Draw();

char text5[100];
sprintf(text5,"#sigma_{J/#psi}= %5.3f #pm %5.3f GeV/c^{2}",functot->GetParameter(6),functot->GetParError(6));
TLatex *l3 = new TLatex(3.5,psifix2->Integral(0.,5.)/binWidth*0.05,text5);
l3->SetTextColor(kBlue);
l3->SetTextSize(0.035);
l3->SetTextFont(42);
l3->Draw();

char text7[100];
sprintf(text7,"#chi^{2}/ndf = %3.2f",chi2/ndf);
TLatex *l5 = new TLatex(4.,psifix2->Integral(0.,5.)/binWidth*0.001,text7);
l5->SetTextColor(kBlue);
l5->SetTextSize(0.038);
l5->SetTextFont(42);
l5->Draw();

TF1 *psipfix = new TF1("psipfix",FuncPsiP,0.,5.,12); 
psipfix->SetParameter(4,functot->GetParameter(4));  
psipfix->SetParameter(5,functot->GetParameter(5));  
psipfix->SetParameter(6,functot->GetParameter(6));  
psipfix->SetParameter(7,functot->GetParameter(7));  
psipfix->SetParameter(8,functot->GetParameter(8)); 
psipfix->SetParameter(9,functot->GetParameter(9)); 
psipfix->SetParameter(10,functot->GetParameter(10)); 
psipfix->SetParameter(11,functot->GetParameter(11)); 
psipfix->SetLineColor(kAzure+6);
psipfix->Draw("same");

TF1 *psipfix2 = new TF1("psipfix2",FuncPsiP2,0.,5.,8); 
psipfix2->SetParameter(0,functot->GetParameter(4));
psipfix2->SetParameter(1,functot->GetParameter(5));
psipfix2->SetParameter(2,functot->GetParameter(6));
psipfix2->SetParameter(3,functot->GetParameter(7)); 
psipfix2->SetParameter(4,functot->GetParameter(8));  
psipfix2->SetParameter(5,functot->GetParameter(9));  
psipfix2->SetParameter(6,functot->GetParameter(10));  
psipfix2->SetParameter(7,functot->GetParameter(11));  
Double_t NPsiP = psipfix->Integral(0.,5.)/binWidth;

 Double_t psipmat[64];
 for(Int_t i=0;i<8;i++){
   for(Int_t j=0;j<8;j++){
     psipmat[8*i+j]=fullmat[52+j+12*i];
   }
 }
 Double_t psippar[8];
 for(int i=0;i<8;i++) {
   psippar[i]=functot->GetParameter(4+i);
 } 
 Double_t ErrPsiPCorrParam = psipfix2->IntegralError(0.,5.,psippar,psipmat)/binWidth;
 Double_t ErrNPsiP=ErrPsiPCorrParam;

 printf("PsiP=%f +/- %f\n",NPsiP,ErrNPsiP);

char text4[100];
sprintf(text4,"N_{#psi}(2S)= %3.0f #pm %2.0f",NPsiP,ErrNPsiP);
TLatex *l2 = new TLatex(3.5,psifix2->Integral(0.,5.)/binWidth*0.006,text4);
l2->SetTextColor(kBlue);
l2->SetTextSize(0.038);
l2->SetTextFont(42);
l2->Draw();

char text40[100];
sprintf(text40,"m_{#psi(2S)}= %5.3f GeV/c^{2}",functot->GetParameter(5)+(3.686-3.097));
TLatex *l20 = new TLatex(3.5,psifix2->Integral(0.,5.)/binWidth*0.004,text40);
l20->SetTextColor(kBlue);
l20->SetTextSize(0.035);
l20->SetTextFont(42);
l20->Draw();

char text50[100];
sprintf(text50,"#sigma_{#psi(2S)}= %5.3f GeV/c^{2}",functot->GetParameter(6)*(functot->GetParameter(5)+(3.686-3.097))/functot->GetParameter(5));
TLatex *l30 = new TLatex(3.5,psifix2->Integral(0.,5.)/binWidth*0.0025,text50);
l30->SetTextColor(kBlue);
l30->SetTextSize(0.035);
l30->SetTextFont(42);
l30->Draw();

 TF1 *bckfix = new TF1("bckfix",FuncBck,FitMin,FitMax,4);  
 bckfix->SetParameter(0,functot->GetParameter(0));  
 bckfix->SetParameter(1,functot->GetParameter(1));  
 bckfix->SetParameter(2,functot->GetParameter(2));  
 bckfix->SetParameter(3,functot->GetParameter(3));  
 bckfix->SetLineColor(kOrange);
 bckfix->Draw("same");
 
 functot->Draw("same");
 
 Double_t Sigma3Min=functot->GetParameter(5)-3.*functot->GetParameter(6);
 Double_t Sigma3Max=functot->GetParameter(5)+3.*functot->GetParameter(6);
 Double_t NPsi3Sigma= psifix->Integral(Sigma3Min,Sigma3Max)/binWidth;
 Double_t NBck3Sigma = bckfix->Integral(Sigma3Min,Sigma3Max)/binWidth;
 Double_t Significance = NPsi3Sigma/TMath::Sqrt(NPsi3Sigma+NBck3Sigma);
 
 printf("S/B (3sigma)=%f\n",NPsi3Sigma/NBck3Sigma);
 printf("S/sqrt(S+B) (3sigma)=%f\n",Significance);

 Double_t Sigma3PsiPMin=(functot->GetParameter(5)+(3.686-3.097))-3.*functot->GetParameter(6)*(functot->GetParameter(5)+(3.686-3.097))/functot->GetParameter(5);
 Double_t Sigma3PsiPMax=(functot->GetParameter(5)+(3.686-3.097))+3.*functot->GetParameter(6)*(functot->GetParameter(5)+(3.686-3.097))/functot->GetParameter(5);
 printf("Sigma3PsiPMin=%f Sigma3PsiPMax=%f\n",Sigma3PsiPMin,Sigma3PsiPMax);
 Double_t NPsiP3Sigma= psipfix->Integral(Sigma3PsiPMin,Sigma3PsiPMax)/binWidth;
 Double_t NBckPsiP3Sigma = bckfix->Integral(Sigma3PsiPMin,Sigma3PsiPMax)/binWidth;
 printf("NPsiP3Sigma=%f ,NBckPsiP3Sigma=%f\n ",NPsiP3Sigma,NBckPsiP3Sigma);
 Double_t SignificancePsiP = NPsiP3Sigma/TMath::Sqrt(NPsiP3Sigma+NBckPsiP3Sigma);
 
 char text60[100];
sprintf(text60,"S/B Psi(2S)(3#sigma) = %3.2f",NPsiP3Sigma/NBckPsiP3Sigma);
TLatex *l40 = new TLatex(3.5,psifix2->Integral(0.,5.)/binWidth*0.2,text60);
l40->SetTextColor(kBlue);
l40->SetTextSize(0.038);
l40->SetTextFont(42);
printf("S/B Psi(2S)(3#sigma) = %5.4f\n",NPsiP3Sigma/NBckPsiP3Sigma);
printf("S/sqrt(S+B) Psi(2S) (3sigma)=%f\n",SignificancePsiP);

 Double_t Significance = NPsi3Sigma/TMath::Sqrt(NPsi3Sigma+NBck3Sigma);

char text6[100];
sprintf(text6,"S/B (3#sigma) = %3.2f",NPsi3Sigma/NBck3Sigma);
TLatex *l4 = new TLatex(3.5,psifix2->Integral(0.,5.)/binWidth*0.15,text6);
l4->SetTextColor(kBlue);
l4->SetTextSize(0.038);
l4->SetTextFont(42);
l4->Draw();

}


Double_t FuncJpsi(Double_t *x, Double_t *par){

  Double_t t = (x[0]-par[5])/par[6];
  if (par[7] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[7]);
  Double_t absAlpha2 = fabs((Double_t)par[9]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[4]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[8]/absAlpha,par[8])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[8]/absAlpha - absAlpha;

    return par[4]*(a/TMath::Power(b - t, par[8]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[10]/absAlpha2,par[10])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[10]/absAlpha2 - absAlpha2;

  return  par[4]*(c/TMath::Power(d + t, par[10]));
  }

  return 0. ;
}

Double_t FuncJpsi2(Double_t *x, Double_t *par){
  /*
par[1] = mean;
par[2] = sigma;
par[3] = alpha;
par[4] = n;
par[0] = Normalization
par[5] = alpha2;
par[6] = n2;
*/

  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[0]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;

    return par[0]*(a/TMath::Power(b - t, par[4]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[6]/absAlpha2 - absAlpha2;

  return  par[0]*(c/TMath::Power(d + t, par[6]));
  }

  return 0. ;
}

Double_t FuncPsiP(Double_t *x, Double_t *par){

  Double_t t = (x[0]-(par[5]+(3.686-3.097)))/(par[6]*(par[5]+(3.686-3.097))/par[5]);
  if (par[7] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[7]);
  Double_t absAlpha2 = fabs((Double_t)par[9]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[4]*par[11]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[8]/absAlpha,par[8])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[8]/absAlpha - absAlpha;

    return par[4]*par[11]*(a/TMath::Power(b - t, par[8]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[10]/absAlpha2,par[10])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[10]/absAlpha2 - absAlpha2;

  return  par[4]*par[11]*(c/TMath::Power(d + t, par[10]));
  }

  return 0. ;
}

Double_t FuncPsiP2(Double_t *x, Double_t *par){

  Double_t t = (x[0]-(par[1]+(3.686-3.097)))/(par[2]*(par[1]+(3.686-3.097))/par[1]);
  if (par[3] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[0]*par[7]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;

    return par[0]*par[7]*(a/TMath::Power(b - t, par[4]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[6]/absAlpha2 - absAlpha2;

  return  par[0]*par[7]*(c/TMath::Power(d + t, par[6]));
  }

  return 0. ;
}

Double_t FuncBck(Double_t *x, Double_t *par){
  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
  Double_t FitBck = par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
  return FitBck;
}


Double_t FuncTot(Double_t *x, Double_t *par){
  return FuncBck(x,par)+FuncJpsi(x,par)+FuncPsiP(x,par);
}