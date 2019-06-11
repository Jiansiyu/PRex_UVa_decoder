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

void simpro(TString ifi0, TString ifi1="", TString ifi2="") {

  Int_t minhit=2; // if hit on x AND on y is less than minhit, the event is skipped
  
  createGradientColor(20);

  TH1F *h1rx = new TH1F("h1rx","single run x profile", 1280, -0.5, 1279.5); // only up to 1024 used
  h1rx->SetLineColor(2);
  h1rx->SetXTitle("strip idx");
  TH1F *h1ry = new TH1F("h1ry","single rune y profile ", 1280, -0.5, 1279.5);
  h1ry->SetLineColor(4);
  h1ry->SetXTitle("strip idx");
  
  TH1F *h1x = new TH1F("h1x","x profile", 1024, -0.5, 1023.5);
  h1x->SetLineColor(2);
  h1x->SetXTitle("x strip");
  TH1F *h1y = new TH1F("h1y","y profile", 1280, -0.5, 1279.5);
  h1y->SetLineColor(4);
  h1y->SetXTitle("y strip");
    
  TH1F *hsx = new TH1F("hsx","x adc distribution vs sample",6,-0.5,5.5);
  hsx->SetLineColor(2);
  hsx->SetXTitle("Sample");
  TH1F *hsy = new TH1F("hsy","y adc distribution vs sample",6,-0.5,5.5);
  hsy->SetLineColor(4);
  hsy->SetXTitle("Sample");
  
  TH1F *hnx = new TH1F("hnx","Hits number distribution on x", 100, -0.5,99.5);
  hnx->SetLineColor(2);
  hnx->SetXTitle("#hits/event");
  TH1F *hny = new TH1F("hny","Hits number distribution on y", 100, -0.5,99.5);
  hny->SetLineColor(4);
  hny->SetXTitle("#hits/event");

  TH1F *hcx = new TH1F("hcx","Single x-Strip Charge", 100, 0, 20000);
  hcx->SetLineColor(2);
  hcx->SetXTitle("Charge");
  TH1F *hcy = new TH1F("hcy","Single y-Strip Charge", 100, 0, 20000);
  hcy->SetLineColor(4);
  hcy->SetXTitle("Charge");
      
  TH2F *h2 = new TH2F("h2","x / y cross corr (Readout side)",1024,-0.5,1023.5,1280,-0.5,1279.5);
  h2->SetXTitle("x-strip");
  h2->SetYTitle("y-strip");
  
  TCanvas *cc = new TCanvas("cc","Summary",1600,900);
  cc->Divide(3,2);
  cc->Update();
  
  TChain *tda = new TChain("GEMHit");
  tda->Add(ifi0.Data());
  tda->Add(ifi1.Data());
  tda->Add(ifi2.Data());

  Int_t currentFile;
  Int_t evtID, nch;
  Int_t strip[2048];
  Int_t planeID[2048];
  Int_t adc[6][2048];
  Int_t detID[2048];

  Int_t ix[2048], iy[2048];
  Int_t nx, ny;

  tda->SetBranchAddress("evtID",&evtID);
  tda->SetBranchAddress("nch",&nch);
  tda->SetBranchAddress("strip",strip);
  tda->SetBranchAddress("planeID",planeID);
  tda->SetBranchAddress("detID",detID);

  for (Int_t i=0;i<6;i++) {
    tda->SetBranchAddress(Form("adc%d",i),&adc[i][0]);
  }

  cc->cd(1);
  tda->Draw("strip");

  Int_t nn = tda->GetEntries();

  Int_t nfiles = tda->GetListOfFiles()->GetEntries();
    
  printf(" Number of files in chain %d\n",nfiles);

  Int_t idxfile=0;
  TString fnameold="";
  
  TCanvas *cde = new TCanvas("cde","Detailed x/y profiles");
  Int_t ncx,ncy;
  ncx = (Int_t) TMath::Sqrt(nfiles);
  ncy = TMath::Ceil(nfiles/ncx);

  cde->Divide(ncx,ncy);
  cde->Update();

  Int_t evtID_old=-1;
  
  for (Int_t i=0;i<nn;i++) {

    tda->GetEntry(i);

    nx = 0;
    ny = 0;

    // associate ch to x or y axis
    
    for (Int_t j=0;j<nch;j++) {
      if (planeID[j]==0) { // x axis
	h1rx->Fill(1023-strip[j]);	
	ix[nx] = j;
	nx++;
      } else {
	h1ry->Fill(strip[j]);	
	iy[ny] = j;
	ny++;
      }
    }

    if ((evtID<evtID_old)||(i==(nn-1))) {
      printf(" Change File %d : %s\n",idxfile,fnameold.Data());
      cde->cd(idxfile+1);
      h1rx->SetTitle(Form("x(red) y(blue) - %s",fnameold.Data()));
      h1rx->DrawCopy("");
      h1ry->DrawCopy("same");
      cde->Update();
      h1rx->Reset();
      h1ry->Reset();
      idxfile++;
    }
    evtID_old = evtID;
    fnameold = tda->GetFile()->GetName();
   
    //    printf("%d findings %d %d\n",i, nx, ny);
    
    // evaluate xcorrelation for matching

    hnx->Fill(nx);
    hny->Fill(ny);
    
    if ((nx<minhit) && (ny<minhit)) { continue; }
    
    for (Int_t jx = 0; jx<nx; jx++) {
      h1x->Fill(1023-strip[ix[jx]]);
      Float_t ax[6];
      Float_t chargex=0;
      for (Int_t jj=0;jj<6;jj++) {
	ax[jj] = adc[jj][ix[jx]];
	hsx->Fill(jj,ax[jj]);
	chargex += ax[jj];
      }
      hcx->Fill(chargex);
      
      for (Int_t jy = 0; jy<ny; jy++) {
	if (jx==0) { h1y->Fill(strip[iy[jy]]); } // avoid double counting
	
	Float_t ay[6];
	Float_t chargey=0;
	for (Int_t jj=0;jj<6;jj++) {
	  ay[jj] = adc[jj][iy[jy]];
	  hsy->Fill(jj,ay[jj]);
	  chargey += ay[jj];
	}
	if (jx==0) { hcy->Fill(chargey); } // avoid double counting
	
	h2->Fill(1023-strip[ix[jx]], strip[iy[jy]], xcor(ax, ay));
	
      }
    }
    
  }

  TText *txt = new TText(0,1100.,"HV");
  cc->cd(1);
  h2->Draw("colz");
  txt->Draw();
  cc->cd(4);
  h1x->Draw();
  cc->cd(2);
  h1y->Draw();
  cc->cd(3);
  hsx->Draw();
  hsy->Draw("same");
  cc->cd(5)->SetLogy();
  hcx->Draw();
  hcy->Draw("same");
  cc->cd(6)->SetLogy();
  hnx->Draw();
  hny->Draw("same");
  cc->Update();

  //  tda->Draw("evtID:Entry$");  
  
}
