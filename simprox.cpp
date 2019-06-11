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
 * Author: E. Cisbani
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
#include <TLine.h>
#include <TFitResult.h>
#include <TProfile2D.h>

// Custom Color Palette
Int_t *fCPalette;
Int_t fNColor;

Int_t createGradientColor(Int_t nb) {

  Double_t Red[3]    = { 1.0, 0.50, 0.00};
  Double_t Green[3]  = { 1.0, 0.50, 0.00};
  Double_t Blue[3]   = { 1.0, 0.25, 0.00};
  Double_t Length[3] = { 0.00, 0.25, 1.00 };
  Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,nb);
  fCPalette = new Int_t[nb];
  fNColor = nb;
  for (int i=0;i<nb;i++) fCPalette[i] = FI+i;

  return 0;

}

/*
 *
 *
 */

Float_t xcor(Float_t x[6], Float_t y[6]) {

  Float_t sx=0;
  Float_t sy=0;

  Float_t mx=0;
  Float_t my=0;

  Float_t cxy=0;
  Float_t rxy=0;
  
  for (Int_t i=0;i<6;i++) {
    mx += x[i];
    my += y[i];
    sx += x[i]*x[i];
    sy += y[i]*y[i];
    cxy += x[i]*y[i];
  }

  // mean (not used)
  mx = mx / 6.;
  my = my / 6.;

  // standard deviation
  sx = TMath::Sqrt((sx - 6*mx*mx)/6.);
  sy = TMath::Sqrt((sy - 6*my*my)/6.);

  // cross correlation coefficient (at 0)
  //  cxy = 0.; // ???
  if ((sx != 0) && (sy != 0)) {
    rxy = (cxy/6. - mx*my) / (sx*sy);
  }

  //  printf(" xcor : %f\n",rxy);
  
  return rxy;
  
};

Float_t mapx(Int_t det, Float_t local) {

  return local;
  
};

Float_t mapy(Int_t det, Float_t local) {

  return ((det % 3)*1280 + local);

};

Float_t mapz(Int_t det) {

  Float_t zdet[12]={759,759,759,0,0,0,1158,1158,1158,404,404,404};
  
  return zdet[det];
  
}

#define NMODULES 12
#define NCHAMBER NMODULES/3

/*
 *  cor_thr : correlation threshold
 *  res_thr : residual threshold (mm)
 *
 */
void simpro(Float_t cor_thr, Float_t res_thr, TString ifi0, TString ifi1="", TString ifi2="") {

  Int_t gmod[12]={9,15,12,13,14,10,18,0,4,16,19,17};

  Float_t x_offset[12]={-1.234,-1.284,-4.639, 9.616, 10.53, 9.807,7.019, 10.89, 14.89, -3.496, -3.241, -2.296};

  Float_t fracxy= 1024./1280./3; // x/y strips ratio

  Float_t effstat = 1000; // define the size of the efficiency histogram bins
  Float_t chi2_thr = 15; // if chi2 linear fit > chi2_thr, track is not considered
  Int_t minhit=1; // if hit on x AND on y is less than minhit, the event is skipped

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


  
  //  createGradientColor(20);

  TH1F *hcx = new TH1F("hcx","Max Correlation x-Hit", NMODULES, -0.5, NMODULES-0.5);
  hcx->SetLineColor(2);
  hcx->SetXTitle("x-Tracking");
  TH1F *hcy = new TH1F("hcy","Max Correlation y-Hit", NMODULES, -0.5, NMODULES-0.5);
  hcy->SetLineColor(4);
  hcy->SetXTitle("y-Tracking");
      
  TH2F *hcor = new TH2F("hcor","Correlation distributions", NMODULES, -0.5, NMODULES-0.5, 100,0.,1.);

  TH1F *hnhit = new TH1F("hnhit","Hit count distribution", NCHAMBER, -0.5, NCHAMBER-0.5);

  TH1F *hchi2x = new TH1F("hchi2x","Chi2/NDF from x tracking fit", 100, 0, 100);
  TH1F *hchi2y = new TH1F("hchi2y","Chi2/NDF from y tracking fit", 100, 0, 100);

  TH1F *hresx[NMODULES]; // difference between predicted (fit) and measured hit position
  TH1F *hresy[NMODULES];

  TH2F *hxy[NMODULES]; // scatter plot x/y of each GEM module
  TProfile2D *hxycor[NMODULES]; // scatter plot x/y of each GEM module, weighter by correlation
  TH2F *hcxy[NCHAMBER]; // scatter plot x/y of each GEM chamber

  TH1F *hxx[NMODULES];
  TH1F *hyy[NMODULES];

  TH2F *hrhit[NCHAMBER]; // number of reconstructed 2D hit from other chambers
  TH2F *hmhit[NCHAMBER]; // number of matching hit between reconstructed and measured, within a given residual distance

  for (Int_t im=0;im<NMODULES;im++) {
    hxy[im] = new TH2F(Form("hxy%d",im),Form("x/y M%2d",gmod[im]), 1024, -0.5,1023.5,1280,-0.5,1279.5);
    //    hxy[im]->SetMarkerStyle(6);
    hxycor[im] = new TProfile2D(Form("hxycor%d",im),Form("x/y M%2d",gmod[im]), 1024/16, -0.5, 1023.5, 1280/16, -0.5, 1280-0.5);
    //    hxycor[im]->SetMarkerStyle(6);
    hxx[im] = new TH1F(Form("hxx%d",im),Form("x 1D hit M%2d",gmod[im]), 1024, -0.5, 1023.5);
    hyy[im] = new TH1F(Form("hyy%d",im),Form("y 1D hit M%2d",gmod[im]), 1280, -0.5, 1280-0.5);
    hresx[im] = new TH1F(Form("hresx%d",im),Form("x residual (mm) M%2d",gmod[im]), 100,-50,50);
    hresy[im] = new TH1F(Form("hresy%d",im),Form("y residual (mm) M%2d",gmod[im]), 100,-50,50);
  }

  for (Int_t ic=0;ic<NCHAMBER;ic++) {
    hcxy[ic] = new TH2F(Form("hcxy%d",ic),Form("x/y Chamber J%d",ic), 1024, -0.5, 1023.5, 1280*3, -0.5, 1280*3-0.5);
    //    hcxy[ic]->SetMarkerStyle(6);
  }
  
  TH1F *hslx = new TH1F("hslx","Slope X",100,-90.,90.);
  TH1F *hsly = new TH1F("hsly","Slope Y",100,-90.,90.);

  TH2F *hnull = new TH2F("hnull","Track x(red), y(blue) vs z (mm)",10,-10,1200,10,-800,800);

  TText *mytxt = new TText();
  
  TCanvas *cc = new TCanvas("cc",Form("Chamber Summary (%s)",ifi0.Data()),1600,900);
  cc->Divide(NCHAMBER,1);
  cc->Update();

  TCanvas *ccp = new TCanvas("ccp",Form("Chamber Profiles (%s)",ifi0.Data()),1600,900);
  ccp->Divide(NCHAMBER,2);
  ccp->Update();

  TCanvas *cef = new TCanvas("cef",Form("Chamber Efficiency (%s)",ifi0.Data()),1600,900);
  cef->Divide(NCHAMBER,1);
  cef->Update();

  TCanvas *cca = new TCanvas("cca","Tracks Slopes");
  cca->Divide(2,2);
  cca->Update();

  /*
  TCanvas *ccb = new TCanvas("ccb","Tracks Samples");
  ccb->Divide(4,4);
  ccb->Update();
  */

  Int_t ncx,ncy;
  ncx = NCHAMBER;
  ncy = 3;
  /*
  TCanvas *cde = new TCanvas("cde","Detailed x/y profiles",1024,1280);
  cde->Divide(ncx,ncy);
  cde->Update();
  */

  TCanvas *ccc = new TCanvas("ccc","Module Correlation Summary",1024,1280);
  ccc->Divide(ncx,ncy);
  ccc->Update();

  TCanvas *cres = new TCanvas("cres","Module x and y residuals",1024,1280);
  cres->Divide(NCHAMBER,6);
  cres->Update();

  TCanvas *cd1dx = new TCanvas("cd1dx","x 1D Hit profiles",1024,1280);
  cd1dx->Divide(ncx,ncy);
  cd1dx->Update();

  TCanvas *cd1dy = new TCanvas("cd1dy","y 1D Hit profiles",1024,1280);
  cd1dy->Divide(ncx,ncy);
  cd1dy->Update();

  printf("NCX NCY %d %d\n",ncx,ncy);
  
  TChain *tda = new TChain("GEMHit");
  tda->Add(ifi0.Data());
  tda->Add(ifi1.Data());
  tda->Add(ifi2.Data());

  Int_t currentFile;
  Int_t evtID, nch;
  Int_t strip[2048];
  Int_t detID[2048];
  Int_t planeID[2048];
  Int_t adc[6][2048];
  
  Int_t ix[NMODULES][2048], iy[NMODULES][2048]; // number of gem modules, strips index
  Int_t nx[NMODULES], ny[NMODULES];

  tda->SetBranchAddress("evtID",&evtID);
  tda->SetBranchAddress("nch",&nch);
  tda->SetBranchAddress("strip",strip);
  tda->SetBranchAddress("planeID",planeID);
  tda->SetBranchAddress("detID",detID);

  for (Int_t i=0;i<6;i++) { // samples
    tda->SetBranchAddress(Form("adc%d",i),&adc[i][0]);
  }

  cc->cd(1);
  tda->Draw("strip");

  Int_t nentries = tda->GetEntries();

  Int_t nfiles = tda->GetListOfFiles()->GetEntries();
    
  printf("Number of files in chain %d, total events %d\n",nfiles, nentries);

  Float_t n2bin = ((Float_t) nentries) / effstat; // max 2D bins to get reasonable statistics
  Int_t nybin = (Int_t) TMath::Sqrt(n2bin/fracxy);
  Int_t nxbin = (Int_t) TMath::Sqrt(n2bin*fracxy);
  Float_t xlim = 0.4 * 1024/2; //mm
  Float_t ylim = 0.4 * 1280*3/2; // mm
  for (Int_t ic=0;ic<NCHAMBER;ic++) {
    hrhit[ic] = new TH2F(Form("hrhit%d",ic),Form("Number reconstructed hits in J%d",ic), nxbin, -xlim, xlim, nybin, -ylim, +ylim); // mm
    hmhit[ic] = new TH2F(Form("hmhit%d",ic),Form("Efficiency: matched/reco hits in J%d",ic), nxbin, -xlim, xlim, nybin, -ylim, +ylim); // mm
  }

  Int_t idxfile=0;
  TString fnameold="";
  
  Int_t evtID_old=-1;
  Int_t itkr=0;
  
  //  nentries=10000;
  for (Int_t i=0;i<nentries;i++) {

    tda->GetEntry(i);

    for (Int_t im=0;im<NMODULES;im++) {
      nx[im] = 0;
      ny[im] = 0;
    }
    
    // associate ch (strip) to x or y axis

    for (Int_t j=0;j<nch;j++) {
      if (nch>2048) {
	printf(" %d %d : out of bound, skip additional channels\n",i,j);
	break;
      }
      Int_t jm = detID[j]; // index of module
      if (planeID[j]==0) { // x axis
	ix[jm][nx[jm]] = j;
	nx[jm]++;
	hxx[jm]->Fill(strip[j]);
      } else {
	iy[jm][ny[jm]] = j;
	ny[jm]++;
	hyy[jm]->Fill(strip[j]);
      }
    }

    // loop on detectors (GEM modules) to get 2D hit by cross-correlation 
    
    Int_t selx[NMODULES], sely[NMODULES]; // GEM modules
    Float_t selc[NMODULES]; // correlation

    Int_t flaghit=0;
    Float_t maxcor=-1;

    for (Int_t im=0;im<NMODULES;im++) {
      selx[im]=-1;
      sely[im]=-1;
      selc[im]=-1;
      maxcor=-1;
      
      if ((nx[im]<minhit) && (ny[im]<minhit)) { continue; }

      for (Int_t jx = 0; jx<nx[im]; jx++) {
	Float_t axa[6];
	Float_t chargex=0;
	for (Int_t jj=0;jj<6;jj++) {
	  axa[jj] = adc[jj][ix[im][jx]];
	  chargex += axa[jj];
	}
      
	for (Int_t jy = 0; jy<ny[im]; jy++) {	
	  Float_t aya[6];
	  Float_t chargey=0;
	  for (Int_t jj=0;jj<6;jj++) {
	    aya[jj] = adc[jj][iy[im][jy]];
	    chargey += aya[jj];
	  }
	
	  Float_t ccor = xcor(axa, aya);
	  
	  if (ccor>maxcor) { // choose candidate 2d-hit as correlation maximum
	    maxcor = ccor;
	    selx[im] = ix[im][jx];
	    sely[im] = iy[im][jy];
	    selc[im] = ccor;
	  }
	  
	}
	
      }

      if (selc[im]>cor_thr) {
      hcor->Fill(im,selc[im]);

      if ((selx[im]>=0) && (sely[im]>=0)) {
	hxy[im]->Fill(strip[selx[im]], strip[sely[im]]);
	hxycor[im]->Fill(strip[selx[im]], strip[sely[im]], selc[im]);
	flaghit=1; // at least 1 hit
      }
      }
      
    } // im

    Float_t cmx[NCHAMBER], cmy[NCHAMBER], cmz[NCHAMBER]; // GEM chamber hits
    Int_t cmidx[NCHAMBER];  // chamber layer index with hit
    Float_t residual[NCHAMBER];
    Int_t midx[NCHAMBER]; // module that fired

    Float_t xx,yy,zz;
    Float_t xyc;
    Int_t count=0;
    
    for (Int_t ic=0;ic<NCHAMBER;ic++) { // chamber
      residual[ic]=1.0e9;
      xx=yy=zz=-1;
      xyc=-1;
      cmidx[ic]=-1;
      for (Int_t iml=0;iml<3;iml++) { // module
	Int_t im = ic*3+iml;
	if (selc[im]>xyc) {
	  xx = mapx(im,strip[selx[im]]);
	  yy = mapy(im,strip[sely[im]]);
	  zz = mapz(im);
	  xyc=selc[im];
	  midx[ic]=im;
	}
      }
      if (zz>=0) {
	if (xyc>cor_thr) { // correlation threshold
	  hcxy[ic]->Fill(xx,yy);
	  cmx[count]=(xx+0.5)*0.4-512*0.4 - x_offset[midx[ic]];
	  cmy[count]=(yy+0.5)*0.4-1280*1.5*0.4;
	  cmz[count]=zz;
	  cmidx[ic]=count;
	  count++;
	}
      }
    }

    /*
     * --------- tracking
     */

    hnhit->Fill(count);

    if (count>2) { // at least hits in 3 different chambers
      for (Int_t h=0;h<NCHAMBER;h++) { // loop on chambers
	Int_t ilc = cmidx[h];
	
	if ((ilc>=0) && (count<4)) { // not enough hits for tracking, at least 3 required
	  continue; // this will not be processed
	}

	TGraph *gx = new TGraph(count,cmz,cmx);
	TGraph *gy = new TGraph(count,cmz,cmy);
	gx->SetMarkerStyle(20);
	gy->SetMarkerStyle(22);

	if (ilc>=0) { // remove chamber from fit
	  gx->RemovePoint(ilc);
	  gy->RemovePoint(ilc);
	}

	TFitResultPtr fitr = gx->Fit("pol1","FQS,EX0"); 
	Double_t chi2x = fitr->Chi2() / fitr->Ndf();
	Double_t ax = fitr->GetParams()[0];
	Double_t bx = fitr->GetParams()[1];
	Double_t slopex = TMath::ATan(bx)*160./3.1415;

	fitr = gy->Fit("pol1","FQS,EX0"); 
	Double_t chi2y = fitr->Chi2() / fitr->Ndf();
	Double_t ay = fitr->GetParams()[0];
	Double_t by = fitr->GetParams()[1];
	Double_t slopey = TMath::ATan(by)*160./3.1415;

	gy->GetFunction("pol1")->SetLineColor(4);

	hchi2x->Fill(chi2x);
	hchi2y->Fill(chi2y);

	if ((chi2x<chi2_thr)&&(chi2y<chi2_thr)) { // accept reasonable tracks

	  Float_t zc = mapz(h*3); // z of first module in chamber is the reference 
	  Float_t xfit = ax + bx * zc; 
	  Float_t yfit = ay + by * zc;

	  hrhit[h]->Fill(xfit, yfit);

	  if (ilc>=0) {
	    Float_t xres = cmx[ilc]-xfit;
	    Float_t yres = cmy[ilc]-yfit;
	    hresx[midx[h]]->Fill(xres);
	    hresy[midx[h]]->Fill(yres);
	    residual[h] = TMath::Sqrt(xres*xres+yres*yres);
	  }

	  if (residual[h]<res_thr) {
	    hmhit[h]->Fill(xfit, yfit);
	  }
	}

	delete gx;
	delete gy;

      } // loop on chambers

	/*
	  if ((chi2x<100)&&(chi2y<100)) {
	  hslx->Fill(slopex);
	  hsly->Fill(slopey);
	  }
	*/

	/*
	  if ((itkr<16)&&(count>=3)&&(chi2x<100)&&(chi2y<100)) {
	  ccb->cd(itkr%16+1);
	  hnull->Draw();
	  gy->Draw("P");
	  gx->Draw("P");
	  mytxt->DrawText(100,600,Form("%.1f %.1f (%d)",chi2x,chi2y,evtID));
	  ccb->Update();
	  itkr++;
	  }
	*/

    }

    if ((i%1000) == 0) { printf("Processing event %d out of %d\n",i, nentries); }

  } // entries
  
  //  cc->cd(1);
  //  hcor->Draw("colz");
  
  for (Int_t ic=0;ic<NCHAMBER;ic++) {
    cc->cd(ic+1);
    hcxy[ic]->Draw("");
    ccp->cd(ic+1);
    hcxy[ic]->ProjectionX()->Draw("");
    ccp->cd(ic+NCHAMBER+1);
    hcxy[ic]->ProjectionY()->Draw("");
    //    cef->cd(ic+1);
    //    hrhit[ic]->Draw("colz");
    cef->cd(ic+1);
    hmhit[ic]->Divide(hrhit[ic]);
    hmhit[ic]->GetZaxis()->SetRangeUser(0,1);
    hmhit[ic]->Draw("colz");
    hrhit[ic]->Draw("same");
  }
  cc->Update();
  ccp->Update();
  cef->Update();

  cca->cd(1);
  hslx->Draw();
  cca->cd(2);
  hsly->Draw();
  cca->cd(3);
  hchi2x->Draw();
  cca->cd(4);
  hchi2y->Draw();
  cca->Update();

  for (Int_t im=0;im<NMODULES;im++) {
    Int_t idx = im/3;
    Int_t idy = 2 - im % 3; // 2, 1, 0
    Int_t idp = idy * NCHAMBER + idx + 1;
    
    //    printf("%d: %d %d %d\n",im, idx,idy,idp);
    //    cde->cd(idp);
    //    hxy[im]->Draw("");
    ccc->cd(idp);
    hxycor[im]->GetZaxis()->SetRangeUser(0,1);
    hxycor[im]->Draw("colz");
    hxy[im]->Draw("same");
    cd1dx->cd(idp);
    hxx[im]->Draw();
    cd1dy->cd(idp);
    hyy[im]->Draw();

    //    gStyle->SetOptStat(1);
    //    gStyle->SetOptFit(1);
    cres->cd(idy*NCHAMBER + idx + 1);
    hresx[im]->Draw();
    hresx[im]->Fit("gaus");
    cres->cd(idy*NCHAMBER + idx + 1 + 3*NCHAMBER);
    hresy[im]->Fit("gaus");
    hresy[im]->Draw();
    cres->Update();
    //    gStyle->SetOptStat(0);
    //    gStyle->SetOptFit(0);
  }

  //  cde->Update();
  ccc->Update();
  cd1dx->Update();
  cd1dy->Update();

  printf(" Tracks detected %d\n",itkr);

}
