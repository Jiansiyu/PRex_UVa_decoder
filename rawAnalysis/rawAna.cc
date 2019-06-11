#include <TTree.h>
#include <TCanvas.h>
#include <stdio.h>
#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>
#include <map>

const int nhits = 4000;
Int_t nch;
Int_t adc[6][nhits];
Int_t det_id[nhits];
Int_t plane_id[nhits];
Int_t strip_id[nhits];

Int_t nfadc;
Int_t fCH[100];
Int_t ftiming[100];
Int_t fadc[100];
Int_t ntdc;
Int_t tCH[100];
Double_t ttiming[100];

//map<int, int> CaloMap;



//0       1        2



//3       5,6      4  
//        7,8


//9,10    13,14   17,18
//11,12   15,16   na,na



void getAdcHitmap()
{
map<int, int> CaloMap;
//fadc 17
  CaloMap[4352]=0; //17-0

  CaloMap[4353]=3;//17-1

  CaloMap[4354]=9;//17-2
  CaloMap[4355]=10;//17-3
  CaloMap[4356]=11;//17-4
  CaloMap[4357]=12;//17-5

  CaloMap[4358]=1;//17-6

  CaloMap[4359]=5;//17-7
  CaloMap[4360]=6;//17-8
  CaloMap[4361]=7;//17-9
  CaloMap[4362]=8;//17-10

  CaloMap[4363]=13;//17-11
  CaloMap[4364]=14;//17-12
  CaloMap[4365]=15;//17-13

  CaloMap[4366]=2;//17-13
  CaloMap[4367]=4;//17-13
    //fadc18
  CaloMap[4608]=16;//18-0

  CaloMap[4609]=17;//18-1
  CaloMap[4610]=18;//18-2

 TCanvas *c_adc = new TCanvas("c_adc","adcDist",10,10,1000,800);
    c_adc->Divide(2,5);

    TCanvas *c_occup = new TCanvas("c_occup","Occupancy",10,10,1000,800);
    c_occup->Divide(2,5);

    TCanvas *c_hit = new TCanvas("c_hit","hitmap",20,20,1000,800);
    c_hit->Divide(2,5);

    TCanvas *c_Sc = new TCanvas("c_Sc","Scintillators",20,20,1000,800);
    c_Sc->Divide(3,6);

    TCanvas *c_Calo_1 = new TCanvas("c_Calo_1","Calorimeters1",20,20,1000,800);
    c_Calo_1->Divide(4,5);
    

    TCanvas *c_Calo_2 = new TCanvas("c_Calo_2","Calorimeters2",20,20,1000,800);
    c_Calo_2->Divide(4,5);
 
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.3);
    gStyle->SetTitleSize(0.09,"t"); 

    TH1F *h_adc[10];
    TH1F *h_occup[10];
    TH1F *h_hitx[5];
    TH1F *h_hity[5];
    //Sc
    TH1F *h_Sc_adc[6];
    TH1F *h_Sc_ftiming[6];
    TH1F *h_Sc_ttiming[6];
    //Calo
    TH1F *h_Calo_adc[21];
    TH1F *h_Calo_ftiming[21];

    //Sc
    for(int i=0;i<6;i++)
    {
      h_Sc_adc[i] = new TH1F(Form("h_Sc_adc_%d",i),Form("h_Sc_adc_%d",i),2000,0,6000);
 
      h_Sc_ftiming[i] = new TH1F(Form("h_Sc_ftiming_%d",i),Form("h_Sc_ftiming_%d",i),50,0,50);

      h_Sc_ttiming[i]= new TH1F(Form("h_Sc_ttiming_%d",i),Form("h_Sc_ttiming_%d",i),200,80,170);

    }

    //Calo
    for(int i=0;i<21;i++)
      {
	h_Calo_adc[i] = new TH1F(Form("h_Calo_adc_%d",i),Form("h_Calo_adc_%d",i),2000,0,20000);
	h_Calo_ftiming[i] = new TH1F(Form("h_Calo_ftiming_%d",i),Form("h_Calo_ftiming_%d",i),50,0,50);
      }

    //GEM
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
    TFile *f = new TFile("/home/danning/hallavme12Decoder/Result/run405_part.root");
    TTree *t = (TTree*)f->Get("GEMHit");

    for(int i=0;i<6;i++)
    {
	t->SetBranchAddress(Form("adc%d", i), &adc[i]);
    }

    t->SetBranchAddress("nch", &nch);
    t->SetBranchAddress("detID", &det_id);
    t->SetBranchAddress("planeID", &plane_id);
    t->SetBranchAddress("strip", &strip_id);

    t->SetBranchAddress("nfadc", &nfadc);
    t->SetBranchAddress("fCH", &fCH);
    t->SetBranchAddress("ftiming", &ftiming);
    t->SetBranchAddress("fadc", &fadc);
    t->SetBranchAddress("ntdc", &ntdc);
    t->SetBranchAddress("tCH", &tCH);
    t->SetBranchAddress("ttiming", &ttiming);

    int N = t->GetEntries();
    cout<<"total entries: "<<N<<endl;
    //N=3000;
    for(int i=0;i<N;i++)
    {
      if(i%100==0) cout<<"Event : "<<i<<endl;
	t->GetEntry(i);
	int p[5][2]={0,0,0,0,0,0,0,0,0,0};

	//GEM
	for(int j=0;j<nch;j++)
	{
	    float aadc = 0.;
	    for(int ii = 0;ii<6;ii++)
	    {
		aadc+=(float)adc[ii][j];
	    }
	    if(aadc<301) continue;
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

	//Scint
	     //fadc

	for(int j=0;j<nfadc;j++)
	  {
	    if((fCH[j]-4608)==5||(fCH[j]-4608)==6)
	      {
	
		h_Sc_adc[fCH[j]-4611]->Fill(fadc[j]);
		
		h_Sc_ftiming[fCH[j]-4611]->Fill(ftiming[j]);
	      }
       	  }
	     //tdc
	for(int j=0;j<ntdc;j++)
	  {
	    
	    if(tCH[j]<=5)
	      {
		h_Sc_ttiming[tCH[j]]->Fill(ttiming[j]);
	      }
	    else if(tCH[j]==15)
	      {
		h_Sc_ttiming[0]->Fill(ttiming[j]);
	      }	      
	  }

	//Caloe
	      //fadc
	for(int j=0;j<nfadc;j++)
	  {
	    if(fCH[j]>4607&&fCH[j]<4611||(fCH[j]>4351&&fCH[j]<4368))
	      {
		h_Calo_adc[CaloMap[fCH[j]]]->Fill(fadc[j]);
		h_Calo_ftiming[CaloMap[fCH[j]]]->Fill(ftiming[j]);
	      }
	    
	  }

    }
    f->Close();



   TFile *resultf = new TFile("test.root","RECREATE");
   //Draw GEM
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
    //Draw Sc
    for(int i=0;i<6;i++)
      {
	c_Sc->cd(i*3+1);
	h_Sc_adc[i]->GetXaxis()->SetLabelSize(0.06);
	h_Sc_adc[i]->GetYaxis()->SetLabelSize(0.06);
	h_Sc_adc[i]->Draw();

	c_Sc->cd(i*3+2);
	h_Sc_ftiming[i]->GetXaxis()->SetLabelSize(0.06);
	h_Sc_ftiming[i]->GetYaxis()->SetLabelSize(0.06);
	h_Sc_ftiming[i]->Draw();

	c_Sc->cd(i*3+3);
	h_Sc_ttiming[i]->GetXaxis()->SetLabelSize(0.06);
	h_Sc_ttiming[i]->GetYaxis()->SetLabelSize(0.06);
	h_Sc_ttiming[i]->Draw();

      }

    //Draw Calo
    for(int i=0;i<5;i++)
      {
	c_Calo_1->cd(1+i*4);
	h_Calo_adc[i]->GetXaxis()->SetLabelSize(0.06);
	h_Calo_adc[i]->GetYaxis()->SetLabelSize(0.06);
	h_Calo_adc[i]->Draw();

	c_Calo_1->cd(3+i*4);
	h_Calo_ftiming[i]->GetXaxis()->SetLabelSize(0.06);
	h_Calo_ftiming[i]->GetYaxis()->SetLabelSize(0.06);
	h_Calo_ftiming[i]->Draw();
      }
    for(int i=5;i<9;i++)
      {
	int j=i-5;
	c_Calo_1->cd(2+j*4);
	h_Calo_adc[i]->GetXaxis()->SetLabelSize(0.06);
	h_Calo_adc[i]->GetYaxis()->SetLabelSize(0.06);
	h_Calo_adc[i]->Draw();

	c_Calo_1->cd(4+j*4);
	h_Calo_ftiming[i]->GetXaxis()->SetLabelSize(0.06);
	h_Calo_ftiming[i]->GetYaxis()->SetLabelSize(0.06);
	h_Calo_ftiming[i]->Draw();
      }



    for(int i=9;i<14;i++)
      {
	int j=i-9;
	c_Calo_2->cd(1+j*4);
	h_Calo_adc[i]->GetXaxis()->SetLabelSize(0.06);
	h_Calo_adc[i]->GetYaxis()->SetLabelSize(0.06);
	h_Calo_adc[i]->Draw();

	c_Calo_2->cd(3+j*4);
	h_Calo_ftiming[i]->GetXaxis()->SetLabelSize(0.06);
	h_Calo_ftiming[i]->GetYaxis()->SetLabelSize(0.06);
	h_Calo_ftiming[i]->Draw();
      }
    for(int i=14;i<19;i++)
      {
	int j=i-14;
	c_Calo_2->cd(2+j*4);
	h_Calo_adc[i]->GetXaxis()->SetLabelSize(0.06);
	h_Calo_adc[i]->GetYaxis()->SetLabelSize(0.06);
	h_Calo_adc[i]->Draw();
	c_Calo_2->cd(4+j*4);
	h_Calo_ftiming[i]->GetXaxis()->SetLabelSize(0.06);
	h_Calo_ftiming[i]->GetYaxis()->SetLabelSize(0.06);
	h_Calo_ftiming[i]->Draw();
      }


    c_Calo_2->Write();
    resultf->Write();

}
