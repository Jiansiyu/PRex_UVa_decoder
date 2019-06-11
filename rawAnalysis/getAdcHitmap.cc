#include <TTree.h>
#include <TCanvas.h>
#include <stdio.h>
#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>

const int nhits = 4000;
Int_t nch;
Int_t adc[6][nhits];
Int_t det_id[nhits];
Int_t plane_id[nhits];
Int_t strip_id[nhits];

void getAdcHitmap()
{
    TCanvas *c_adc = new TCanvas("c_adc","adcDist",10,10,1000,800);
    c_adc->Divide(2,5);

    TCanvas *c_occup = new TCanvas("c_occup","Occupancy",10,10,1000,800);
    c_occup->Divide(2,5);

    TCanvas *c_hit = new TCanvas("c_hit","hitmap",20,20,1000,800);
    c_hit->Divide(2,5);

 
    TH1F *h_adc[10];
    TH1F *h_occup[10];
    TH1F *h_hitx[5];
    TH1F *h_hity[5];

    for(int i=1;i<6;i++)
    {
      h_hitx[i-1]= new TH1F(Form("h_hitx_%d",(i-1)),Form("h_hitx_%d_%d",(i-1)),768,0,1536);
      h_hity[i-1]= new TH1F(Form("h_hity_%d",(i-1)),Form("h_hity_%d_%d",(i-1)),640,0,1280);
      for(int iplane=0;iplane<2;iplane++)
	{
	  h_adc[2*i-2+iplane]= new TH1F(Form("h_adc_%d_%d",(i-1),iplane),Form("h_adc_%d_%d",(i-1),iplane),1000,0,10000);
	  h_occup[2*i-2+iplane]= new TH1F(Form("h_occup_%d_%d",(i-1),iplane),Form("h_occup_%d_%d",(i-1),iplane),100,0,100);
       	}
    }


    //gStyle->SetOptFit(1);
    TFile *f = new TFile("../Result/test.root");
    TTree *t = (TTree*)f->Get("GEMHit");

    for(int i=0;i<6;i++)
    {
	t->SetBranchAddress(Form("adc%d", i), &adc[i]);
    }

    t->SetBranchAddress("nch", &nch);
    t->SetBranchAddress("detID", &det_id);
    t->SetBranchAddress("planeID", &plane_id);
    t->SetBranchAddress("strip", &strip_id);

    int N = t->GetEntries();
    //N=60000;
    cout<<"total entries: "<<N<<endl;

    for(int i=0;i<N;i++)
    {
      if(i%100==0) cout<<"Event : "<<i<<endl;
	t->GetEntry(i);
	int p[5][2]={0,0,0,0,0,0,0,0,0,0};
	for(int j=0;j<nch;j++)
	{
	    float aadc = 0.;
	    for(int ii = 0;ii<6;ii++)
	    {
	      if(aadc<adc[ii][j])
		{
		  aadc=adc[ii][j];
		}
		//	aadc+=(float)adc[ii][j];
	    }
	    h_adc[2*det_id[j]+plane_id[j]] -> Fill(aadc);
	    if(plane_id[j]==0)h_hitx[det_id[j]] -> Fill(strip_id[j]);
	    if(plane_id[j]==1)h_hity[det_id[j]] -> Fill(strip_id[j]);
	    p[det_id[j]][plane_id[j]]+=1;
	}

	for(int ii=1;ii<6;ii++)
	  {
	    for(int iplane=0;iplane<2;iplane++)
	      {
		h_occup[2*ii-2+iplane] -> Fill(p[ii-1][iplane]);
	      }
	  }
    }
    f->Close();



   TFile *resultf = new TFile("raw316_50kTest.root","RECREATE");

    for(int i=1;i<6;i++)
    {
	for(int iplane=0;iplane<2;iplane++)
	{
	    c_adc->cd(2*i-1+iplane);
	    h_adc[2*i-2+iplane]->Draw();
	    h_adc[2*i-2+iplane]->Write();

	    c_occup->cd(2*i-1+iplane);
	    h_occup[2*i-2+iplane]->Draw();
	    h_occup[2*i-2+iplane]->Write();

	    c_hit->cd(2*i-1+iplane);
	    if(iplane==0){h_hitx[i-1]->Draw();h_hitx[i-1]->Write();}
	    if(iplane==1){h_hity[i-1]->Draw();h_hity[i-1]->Write();}
	}
    }
    resultf->Write();

}
