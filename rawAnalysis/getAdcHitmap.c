#include <TCanvas.h>
#include <stdio.h>
#include <Riostream.h>

void getAdcHitmap()
{

 TH1::AddDirectory(0);
 //gStyle->SetOptFit(1);

  TH1F *h_adc_clone[10];
  TH1F *h_hit_clone[10];

 
  TCanvas *c_adc = new TCanvas("c_occup","Occupancy",10,10,1000,800);
  c_adc->Divide(2,5);
 
  
  TCanvas *c_hit = new TCanvas("c_hit","hitmap",20,20,1000,800);
  c_hit->Divide(2,5);

  TFile *f = new TFile("../Result/test.root");
  TTree *t = (TTree*)f->Get("GEMHit");


  TH1F *h_adc[10];
  TH1F *h_hit[10];



  for(int i=1;i<6;i++)
    {
      for(int iplane=0;iplane<2;iplane++)
	{
	  h_adc[2*i-2+iplane]= new TH1F(Form("h_adc_%d_%d",(i-1),iplane),Form("h_adc_%d_%d",(i-1),iplane),1000,0,10000);
	  h_hitx[i-1]= new TH1F(Form("h_hitx_%d",(i-1)),Form("h_hitx_%d_%d",(i-1)),768,0,1536);
	  h_hity[i-1]= new TH1F(Form("h_hity_%d",(i-1)),Form("h_hity_%d_%d",(i-1)),640,0,1280);

	  TCut *cut = new TCut(Form("detID==%d&&planeID==%d",(i-1),iplane));

	  c_adc->cd(2*i-1+iplane);
	  gStyle->SetStatW(0.3);
	  gStyle->SetStatH(0.3);
	  //GEMHit->Draw(Form("(adc0+adc2+adc3+adc4+adc5)>>h_adc_%d_%d(1000,0,10000)",(i-1),iplane),*cut);
	   GEMHit->Project(Form("h_adc_%d_%d",(i-1),iplane),"adc0+adc2+adc3+adc4+adc5",*cut);
	   //h_adc[2*i-2+iplane]->Draw();
	   cout<<h_adc[2*i-2+iplane]->GetEntries()<<endl;
	   h_adc_clone[2*i-2+iplane] = (TH1F*)h_adc[2*i-2+iplane]->Clone(Form("h_adc_clone_%d_%d",(i-1),iplane));
	   cout<<"clone: "<<h_adc_clone[2*i-2+iplane]->GetEntries()<<endl;
	  //  h_adc[2*i-2+iplane]->Fit("landau","","",20,4000);

	  c_hit->cd(2*i-1+iplane);
      	  GEMHit->Project(Form("h_hit_%d_%d",(i-1),iplane),"strip",*cut);
	  h_hit[2*i-2+iplane]->Draw();
	  //	  GEMHit->Draw(Form("(adc0+adc2+adc3+adc4+adc5)>>h_adc_%d_%d(1000,0,10000)",(i-1),iplane),*cut);

	}
    }
  

   c_adc->Update();
  c_hit->Update();
  getchar();




  cut->Delete();
  //c_adc->Delete();
  f->Close();

  
  for(int i=1;i<6;i++)
    {
      for(int iplane=0;iplane<2;iplane++)
	{
	  c_adc->cd(2*i-1+iplane);
	  h_adc_clone[2*i-2+iplane]->Draw();
	}
    }
  c_adc->Update();
  getchar();
  


}
