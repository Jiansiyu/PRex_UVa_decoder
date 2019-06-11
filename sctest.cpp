/*
 * Attempt to analyse the single module test root data
 *
 * the root file must be processed by mpd*_decoder or similar, with pedestal subtraction
 * 
 *
 * run root
 * .L simpro.cpp+
 *
 * simpro("test_1.root")
 *
 * where test_1.root is the output of mdp*_decoder ..
 *
 * the filename can contain wildcard to load multiple files in one shot
 * two additional filename parameters can be provided to simpro,
 *  to add root file to the analysis
 *
 */


#include <stdio.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCut.h>
#include <TMath.h>
#include <TColor.h>
#include <TChain.h>
#include <TLegend.h>
#include <TText.h>

// Custom Color Palette
Int_t *fCPalette;
Int_t fNColor;

Int_t createGradientColor(Int_t nb) {

  Double_t Red[3]    = { 1.0, 0.50, 0.00};
  Double_t Green[3]  = { 1.0, 0.50, 0.00};
  Double_t Blue[3]   = { 1.0, 0.25, 0.00};
  Double_t Length[3] = { 0.00, 0.50, 1.00 };
  Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,nb);
  fCPalette = new Int_t[nb];
  fNColor = nb;
  for (int i=0;i<nb;i++) fCPalette[i] = FI+i;

  return 0;

}


Float_t xcor(Float_t x[6], Float_t y[6]) {

  Float_t sx=0;
  Float_t sy=0;

  Float_t mx=0;
  Float_t my=0;

  Float_t cxy=0;
  
  for (Int_t i=0;i<6;i++) {
    mx += x[i];
    my += y[i];
    sx += x[i]*x[i];
    sy += y[i]*y[i];
    cxy += x[i]*y[i];
  }

  // standard deviation
  sx = TMath::Sqrt((sx - mx)/5.);
  sy = TMath::Sqrt((sy - my)/5.);

  // cross correlation coefficient (at 0)
  cxy = 0.;
  if ((sx != 0) && (sy != 0)) {
    cxy = (cxy - mx*my/6.)/6./(sx*sy);
  }

  // mean (not used)
  mx = mx / 6.;
  my = my / 6.;

  //  printf(" xcor : %f\n",cxy);
  
  return -cxy;
  
}

void sctest(TString ifi0) {

  TCanvas *c1 = new TCanvas("c1",Form("Pedestal RMS %s",ifi0.Data()));
  c1->Divide(4,4);
  c1->Update();

  TH1F *hmean = new TH1F("hmean","Average RMS of each APV",16,-0.5,15.5);

  TFile *f1 = new TFile(ifi0);
  f1->ls();

  TH1F *h1=NULL;
  for (Int_t i=0;i<16;i++) {
    h1 = (TH1F*) f1->Get(Form("PedestalRMS_mpd_1_ch_%d",i));
    c1->cd(i+1);
    if (h1 != NULL) {
      h1->DrawCopy();
      Float_t mean=0;
      for (Int_t j=0;j<128;j++) {
	mean += h1->GetBinContent(j+1);
      }
      printf(" Ch %d : RMS mean = %f\n",i, mean/128.);
      hmean->Fill(i,mean/128.);
    }
    h1 = NULL;
  }

  c1->cd(16);
  hmean->Draw();

  f1->Close();
}
