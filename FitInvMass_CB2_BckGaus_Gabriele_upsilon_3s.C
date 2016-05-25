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

Double_t FuncUpsilon(Double_t *, Double_t *);
Double_t FuncUpsilon2(Double_t *, Double_t *);
Double_t FuncUpsilonP(Double_t *, Double_t *);
Double_t FuncUpsilonP2(Double_t *, Double_t *);
Double_t FuncBck(Double_t *, Double_t *);
Double_t FuncTot(Double_t *, Double_t *);

Double_t NormUpsilonP=0;

Double_t UpsilonMass=9.460;
Double_t Upsilon2SMass=10.023;
Double_t Upsilon3SMass=10.355;

void FitInvMass_CB2_BckGaus_Gabriele_upsilon_3s(TString inputFile, TString histoName, Int_t rebinFactor=1, Double_t FitMin=7., Double_t FitMax=12., char *period = "LHC15o", char *tailsfix="yes"){

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
// choose the JUpsilon tails
//-----------------------------------------------------------------------

Double_t cbalpha1; 
Double_t cbn1;
Double_t cbalpha2;
Double_t cbn2;

cbalpha1=0.93; 
cbn1=2.27;
cbalpha2=2.09;
cbn2=2.78;  

printf("\nDouble CB tails:\n");   
printf("\nAlpha1=%f N1=%f Alpha2=%f N2=%f\n",cbalpha1,cbn1,cbalpha2,cbn2);

//-----------------------------------------------------------------------
// prepare plot
//-----------------------------------------------------------------------
histo->Rebin(rebinFactor);

histo->GetXaxis()->SetTitle("M_{#mu#mu} (GeV/c^{2})");
histo->GetYaxis()->SetTitle("dN/dM_{#mu#mu}");
Double_t binWidth= histo->GetBinWidth(1);

Double_t minUpsilonMass=UpsilonMass-0.054*4;
Double_t maxUpsilonMass=Upsilon2SMass+0.054*4;

TCanvas *c = new TCanvas("c","c",20,20,600,600);
c->cd(1);
// set to 0 the entries of the JUpsilon region (minUpsilonMass-maxUpsilonMass)
TH1D *histoForBck = (TH1D*) histo->Clone("histoForBck");
Double_t xxmin=minUpsilonMass/histo->GetBinWidth(1);
Double_t xxmax=maxUpsilonMass/histo->GetBinWidth(1);
histoForBck->GetYaxis()->SetRangeUser(0.00001,2000.);

for(int i=xxmin;i<xxmax;i++){
  histoForBck->SetBinContent(i+1,0);
  histoForBck->SetBinError(i+1,0);
}
//---------------------------------------------------------------------
// 1st step: Bck fit
//---------------------------------------------------------------------
//Double_t par[12]={90000.,0.46,1.14,0.06,5000.,UpsilonMass,7.0e-02,1.00011e+00,Upsilon2SMass,1.68359e+00,Upsilon2SMass,0.01}; 
Double_t par[13]={3.e+06,-1.,1.4e+03,-0.2,98.,UpsilonMass,1.6e-01,1.00011e+00,3.70078e+00,1.68359e+00,3.63002e+00,0.01,0.01}; 
TF1 *funcbck;
ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(300000);
funcbck = new TF1("funcbck",FuncBck,FitMin,FitMax,4); 
funcbck->SetParameter(0,par[0]);
funcbck->SetParameter(1,par[1]);
funcbck->SetParLimits(1,-10.,0.);
funcbck->SetParameter(2,par[2]);
funcbck->SetParameter(3,par[3]);
funcbck->SetParLimits(3,-10.,0.);
funcbck->SetLineColor(kOrange);
funcbck->SetLineWidth(2);
printf("\n--- Step 1: Fitting the bck1+bck2 ---\n");
histoForBck->Fit(funcbck,"RM");

//---------------------------------------------------------------------
// 2nd step: JUpsilon fit
//---------------------------------------------------------------------

TF1 *funcUpsilon = new TF1("funcUpsilon",FuncUpsilon,minUpsilonMass,maxUpsilonMass,12);
//funcUpsilon->SetParameters(par);
funcUpsilon->FixParameter(0,funcbck->GetParameter(0));
funcUpsilon->FixParameter(1,funcbck->GetParameter(1));
funcUpsilon->FixParameter(2,funcbck->GetParameter(2));
funcUpsilon->FixParameter(3,funcbck->GetParameter(3));
funcUpsilon->SetParameter(5,par[5]);
funcUpsilon->SetParLimits(5,9.,10.);
funcUpsilon->FixParameter(6,par[6]);
funcUpsilon->FixParameter(7,cbalpha1);
funcUpsilon->FixParameter(8,cbn1);
funcUpsilon->FixParameter(9,cbalpha2);
funcUpsilon->FixParameter(10,cbn2);
funcUpsilon->FixParameter(11,par[11]);
funcUpsilon->SetParLimits(6,0.,1.);  //fixing sigma defined positive!

funcUpsilon->SetLineColor(kGreen);
funcUpsilon->SetLineWidth(2);
printf("\n--- Step 2: Fitting the Upsilon ---\n");
histo->Fit(funcUpsilon,"ELR");

//---------------------------------------------------------------------
// 3rd step: total fit
//---------------------------------------------------------------------
TF1 *functot = new TF1("functot",FuncTot,FitMin,FitMax,13); 
functot->SetParameter(0,funcbck->GetParameter(0));   
functot->SetParameter(1,funcbck->GetParameter(1));   
functot->SetParameter(2,funcbck->GetParameter(2));   
functot->SetParameter(3,funcbck->GetParameter(3));   
functot->SetParameter(4,funcUpsilon->GetParameter(4));   
functot->SetParameter(5,UpsilonMass); 
functot->SetParameter(6,0.1);

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
functot->SetParameter(12,0.01);   
functot->SetParLimits(12,1.e-04,.1);   

functot->SetLineColor(kBlue);
functot->SetLineWidth(2);

printf("\n--- Step 3: Fitting total spectrum ---\n");
TFitResultPtr r = histo->Fit(functot,"RLS");

TMatrixDSym cov = r->GetCovarianceMatrix();
cov.Print();

TMatrixDSym cov2 = r->GetCorrelationMatrix();
cov2.Print();

Double_t *fullmat;
fullmat = cov.GetMatrixArray();
Double_t Upsilonmat[49];
for(Int_t i=0;i<7;i++){
  for(Int_t j=0;j<7;j++){
    Upsilonmat[7*i+j]=fullmat[56+j+13*i];
  }
}

histo->GetXaxis()->SetRangeUser(FitMin,FitMax);
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


TPaveText *paveText=new TPaveText(.65,.6,0.9,.9,"NDC");
paveText->SetTextColor(kBlue);
paveText->SetTextSize(0.020);
paveText->SetTextFont(42);
paveText->SetShadowColor(kWhite);

//------------------------------------------------------
// plot fitting functions and compute uncertainties
//------------------------------------------------------
TF1 *Upsilonfix = new TF1("Upsilonfix",FuncUpsilon,FitMin,FitMax,12); 
Upsilonfix->SetParameter(0,functot->GetParameter(0));
Upsilonfix->SetParameter(1,functot->GetParameter(1));
Upsilonfix->SetParameter(2,functot->GetParameter(2));
Upsilonfix->SetParameter(3,functot->GetParameter(3)); 
Upsilonfix->SetParameter(4,functot->GetParameter(4));  
Upsilonfix->SetParameter(5,functot->GetParameter(5));  
Upsilonfix->SetParameter(6,functot->GetParameter(6));  
Upsilonfix->SetParameter(7,functot->GetParameter(7));  
Upsilonfix->SetParameter(8,functot->GetParameter(8));  
Upsilonfix->SetParameter(9,functot->GetParameter(9));  
Upsilonfix->SetParameter(10,functot->GetParameter(10));  
Upsilonfix->SetParameter(11,functot->GetParameter(11));  
Upsilonfix->SetLineColor(kGreen);
Upsilonfix->Draw("same");
Double_t NUpsilon=Upsilonfix->Integral(FitMin,FitMax)/binWidth;

TF1 *Upsilonfix2 = new TF1("Upsilonfix2",FuncUpsilon2,FitMin,FitMax,7); 
Upsilonfix2->SetParameter(0,functot->GetParameter(4));
Upsilonfix2->SetParameter(1,functot->GetParameter(5));
Upsilonfix2->SetParameter(2,functot->GetParameter(6));
Upsilonfix2->SetParameter(3,functot->GetParameter(7)); 
Upsilonfix2->SetParameter(4,functot->GetParameter(8));  
Upsilonfix2->SetParameter(5,functot->GetParameter(9));  
Upsilonfix2->SetParameter(6,functot->GetParameter(10));  
Double_t Upsilonpar[7];
for(int i=0;i<7;i++) {
  Upsilonpar[i]=functot->GetParameter(4+i);
} 
Double_t ErrUpsilonCorrParam = Upsilonfix2->IntegralError(FitMin,FitMax,Upsilonpar,Upsilonmat)/binWidth;
Double_t ErrNUpsilon=ErrUpsilonCorrParam;
printf("\nUpsilon=%f \n",Upsilonfix2->Integral(FitMin,FitMax)/binWidth);
printf("ErrUpsilon=%4.3f (computed using corr. params) %f, %f\n",ErrUpsilonCorrParam, Upsilonfix2->IntegralError(FitMin,FitMax,Upsilonpar,Upsilonmat),binWidth);

paveText->AddText(Form("N_{#Upsilon}= %3.0f #pm %2.0f",NUpsilon,ErrUpsilonCorrParam));

paveText->AddText(Form("m_{#Upsilon}= %5.3f #pm %5.3f GeV/c^{2}",functot->GetParameter(5),functot->GetParError(5)));

paveText->AddText(Form("#sigma_{#Upsilon}= %5.3f #pm %5.3f GeV/c^{2}",functot->GetParameter(6),functot->GetParError(6)));

paveText->AddText(Form("#chi^{2}/ndf = %3.2f",chi2/ndf));

TF1 *Upsilonpfix = new TF1("Upsilonpfix",FuncUpsilonP,FitMin,FitMax,12); 
Upsilonpfix->SetParameter(4,functot->GetParameter(4));  
Upsilonpfix->SetParameter(5,functot->GetParameter(5));  
Upsilonpfix->SetParameter(6,functot->GetParameter(6));  
Upsilonpfix->SetParameter(7,functot->GetParameter(7));  
Upsilonpfix->SetParameter(8,functot->GetParameter(8)); 
Upsilonpfix->SetParameter(9,functot->GetParameter(9)); 
Upsilonpfix->SetParameter(10,functot->GetParameter(10)); 
Upsilonpfix->SetParameter(11,functot->GetParameter(11)); 
Upsilonpfix->SetLineColor(kAzure+6);
Upsilonpfix->Draw("same");

TF1 *Upsilonpfix2 = new TF1("Upsilonpfix2",FuncUpsilonP2,FitMin,FitMax,8); 
Upsilonpfix2->SetParameter(0,functot->GetParameter(4));
Upsilonpfix2->SetParameter(1,functot->GetParameter(5));
Upsilonpfix2->SetParameter(2,functot->GetParameter(6));
Upsilonpfix2->SetParameter(3,functot->GetParameter(7)); 
Upsilonpfix2->SetParameter(4,functot->GetParameter(8));  
Upsilonpfix2->SetParameter(5,functot->GetParameter(9));  
Upsilonpfix2->SetParameter(6,functot->GetParameter(10));  
Upsilonpfix2->SetParameter(7,functot->GetParameter(11));  
Double_t NUpsilonP = Upsilonpfix->Integral(FitMin,FitMax)/binWidth;

 Double_t Upsilonpmat[64];
 for(Int_t i=0;i<8;i++){
   for(Int_t j=0;j<8;j++){
     Upsilonpmat[8*i+j]=fullmat[56+j+13*i];
   }
 }
 Double_t Upsilonppar[8];
 for(int i=0;i<8;i++) {
   Upsilonppar[i]=functot->GetParameter(4+i);
 } 
 Double_t ErrUpsilonPCorrParam = Upsilonpfix2->IntegralError(FitMin,FitMax,Upsilonppar,Upsilonpmat)/binWidth;
 Double_t ErrNUpsilonP=ErrUpsilonPCorrParam;

 printf("UpsilonP=%f +/- %f\n",NUpsilonP,ErrNUpsilonP);

paveText->AddText(Form("N_{#Upsilon(2S)}= %3.0f #pm %2.0f",NUpsilonP,ErrNUpsilonP));

paveText->AddText(Form("m_{#Upsilon(2S)}= %5.3f GeV/c^{2}",functot->GetParameter(5)+(Upsilon2SMass-UpsilonMass)));

paveText->AddText(Form("#sigma_{#Upsilon(2S)}= %5.3f GeV/c^{2}",functot->GetParameter(6)*(functot->GetParameter(5)+(Upsilon2SMass-UpsilonMass))/functot->GetParameter(5)));

TF1 *Upsilonpfix3 = new TF1("Upsilonpfix3",FuncUpsilon3S,FitMin,FitMax,13); 
Upsilonpfix3->SetParameter(4,functot->GetParameter(4));  
Upsilonpfix3->SetParameter(5,functot->GetParameter(5));  
Upsilonpfix3->SetParameter(6,functot->GetParameter(6));  
Upsilonpfix3->SetParameter(7,functot->GetParameter(7));  
Upsilonpfix3->SetParameter(8,functot->GetParameter(8)); 
Upsilonpfix3->SetParameter(9,functot->GetParameter(9)); 
Upsilonpfix3->SetParameter(10,functot->GetParameter(10)); 
Upsilonpfix3->SetParameter(12,functot->GetParameter(12)); 
Upsilonpfix3->SetLineColor(kRed);
Upsilonpfix3->Draw("same");

TF1 *Upsilonfix23 = new TF1("Upsilonfix23",FuncUpsilon3S2,FitMin,FitMax,8); 
Upsilonfix23->SetParameter(0,functot->GetParameter(4));
Upsilonfix23->SetParameter(1,functot->GetParameter(5));
Upsilonfix23->SetParameter(2,functot->GetParameter(6));
Upsilonfix23->SetParameter(3,functot->GetParameter(7)); 
Upsilonfix23->SetParameter(4,functot->GetParameter(8));  
Upsilonfix23->SetParameter(5,functot->GetParameter(9));  
Upsilonfix23->SetParameter(6,functot->GetParameter(10));  
Upsilonfix23->SetParameter(7,functot->GetParameter(12));  
Double_t NUpsilon3S = Upsilonpfix3->Integral(FitMin,FitMax)/binWidth;

 Double_t Upsilon3Smat[64];
 for(Int_t i=0;i<8;i++){
   for(Int_t j=0;j<8;j++){
     Int_t i2=i;
     Int_t j2=j;
     if(j==7) j2=8;
     if(i==7) i2=8;
     Upsilon3Smat[8*i+j]=fullmat[56+j2+13*i2];
   }
 }

 Double_t Upsilon3Spar[8];
 for(int i=0;i<8;i++) {
  Int_t i2=i;
  if(i==7) i2=8;
   Upsilon3Spar[i]=functot->GetParameter(4+i2);
 } 

 Double_t ErrUpsilon3SCorrParam = Upsilonfix23->IntegralError(FitMin,FitMax,Upsilon3Spar,Upsilon3Smat)/binWidth;
 Double_t ErrNUpsilon3S=ErrUpsilon3SCorrParam;

printf("Upsilon3S=%f +/- %f\n",NUpsilon3S,ErrNUpsilon3S);

paveText->AddText(Form("N_{#Upsilon(3S)}= %5.4f #pm %2.0f",NUpsilon3S,ErrNUpsilon3S));

paveText->AddText(Form("m_{#Upsilon(3S)}= %5.3f GeV/c^{2}",functot->GetParameter(5)+(Upsilon3SMass-UpsilonMass)));

paveText->AddText(Form("#sigma_{#Upsilon(3S)}= %5.3f GeV/c^{2}",functot->GetParameter(6)*(functot->GetParameter(5)+(Upsilon3SMass-UpsilonMass))/functot->GetParameter(5)));

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
 Double_t NUpsilon3Sigma= Upsilonfix->Integral(Sigma3Min,Sigma3Max)/binWidth;
 Double_t NBck3Sigma = bckfix->Integral(Sigma3Min,Sigma3Max)/binWidth;
 Double_t Significance = NUpsilon3Sigma/TMath::Sqrt(NUpsilon3Sigma+NBck3Sigma);
 
 printf("S/B (3sigma)=%f\n",NUpsilon3Sigma/NBck3Sigma);
 printf("S/sqrt(S+B) (3sigma)=%f\n",Significance);

 Double_t Sigma3UpsilonPMin=(functot->GetParameter(5)+(Upsilon2SMass-UpsilonMass))-3.*functot->GetParameter(6)*(functot->GetParameter(5)+(Upsilon2SMass-UpsilonMass))/functot->GetParameter(5);
 Double_t Sigma3UpsilonPMax=(functot->GetParameter(5)+(Upsilon2SMass-UpsilonMass))+3.*functot->GetParameter(6)*(functot->GetParameter(5)+(Upsilon2SMass-UpsilonMass))/functot->GetParameter(5);
 printf("Sigma3UpsilonPMin=%f Sigma3UpsilonPMax=%f\n",Sigma3UpsilonPMin,Sigma3UpsilonPMax);
 Double_t NUpsilonP3Sigma= Upsilonpfix->Integral(Sigma3UpsilonPMin,Sigma3UpsilonPMax)/binWidth;
 Double_t NBckUpsilonP3Sigma = bckfix->Integral(Sigma3UpsilonPMin,Sigma3UpsilonPMax)/binWidth;
 printf("NUpsilonP3Sigma=%f ,NBckUpsilonP3Sigma=%f\n ",NUpsilonP3Sigma,NBckUpsilonP3Sigma);
 Double_t SignificanceUpsilonP = NUpsilonP3Sigma/TMath::Sqrt(NUpsilonP3Sigma+NBckUpsilonP3Sigma);

printf("S/B Upsilon(2S)(3#sigma) = %5.4f\n",NUpsilonP3Sigma/NBckUpsilonP3Sigma);
printf("S/sqrt(S+B) Upsilon(2S) (3sigma)=%f\n",SignificanceUpsilonP);

 Double_t Significance = NUpsilon3Sigma/TMath::Sqrt(NUpsilon3Sigma+NBck3Sigma);

paveText->AddText(Form("S/B (3#sigma) = %3.2f",NUpsilon3Sigma/NBck3Sigma));

histo->GetYaxis()->SetRangeUser(0.001,12000.);
paveText->Draw();

gPad->Update();
gPad->Modified();

}


Double_t FuncUpsilon(Double_t *x, Double_t *par){

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

Double_t FuncUpsilon2(Double_t *x, Double_t *par){
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

Double_t FuncUpsilonP(Double_t *x, Double_t *par){

  Double_t t = (x[0]-(par[5]+(Upsilon2SMass-UpsilonMass)))/(par[6]*(par[5]+(Upsilon2SMass-UpsilonMass))/par[5]);
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

Double_t FuncUpsilonP2(Double_t *x, Double_t *par){

  Double_t t = (x[0]-(par[1]+(Upsilon2SMass-UpsilonMass)))/(par[2]*(par[1]+(Upsilon2SMass-UpsilonMass))/par[1]);
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

Double_t FuncUpsilon3S(Double_t *x, Double_t *par){

  Double_t t = (x[0]-(par[5]+(Upsilon3SMass-UpsilonMass)))/(par[6]*(par[5]+(Upsilon3SMass-UpsilonMass))/par[5]);
  if (par[7] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[7]);
  Double_t absAlpha2 = fabs((Double_t)par[9]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core

  {
    return par[4]*par[12]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail

  {
    Double_t a =  TMath::Power(par[8]/absAlpha,par[8])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[8]/absAlpha - absAlpha;

    return par[4]*par[12]*(a/TMath::Power(b - t, par[8]));
  }

  if (t >= absAlpha2) //right tail

  {

   Double_t c =  TMath::Power(par[10]/absAlpha2,par[10])*exp(-0.5*absAlpha2*absAlpha2);
   Double_t d = par[10]/absAlpha2 - absAlpha2;

  return  par[4]*par[12]*(c/TMath::Power(d + t, par[10]));
  }

  return 0. ;
}

Double_t FuncUpsilon3S2(Double_t *x, Double_t *par){

  Double_t t = (x[0]-(par[1]+(Upsilon3SMass-UpsilonMass)))/(par[2]*(par[1]+(Upsilon3SMass-UpsilonMass))/par[1]);
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
  //Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
  //Double_t FitBck = par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
  Double_t FitBck = par[0]*TMath::Exp(x[0]*par[1])+par[2]*TMath::Exp(x[0]*par[3]);
  return FitBck;
}


Double_t FuncTot(Double_t *x, Double_t *par){
  return FuncBck(x,par)+FuncUpsilon(x,par)+FuncUpsilonP(x,par)+FuncUpsilon3S(x,par);
}