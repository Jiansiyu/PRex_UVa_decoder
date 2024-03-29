#include <arpa/inet.h>
#include <assert.h>
#include <utility>

#include "raw_decoder.h"

#include <algorithm> 
#include <fstream>
#include <iostream>
#include <utility>

#include <TH1F.h>
#include <TCanvas.h>
using namespace std;

//==========================================================================
RawDecoder::RawDecoder(unsigned int * buffer, int n)
{ 
  mAPVRawSingleEvent.clear();
  mAPVRawHisto.clear();

  fBuf = n;
  buf = new unsigned int[fBuf];
  for(int i=0;i<fBuf;i++)
  {
    buf[i] = buffer[i];
  }
  Decode();
}

//==========================================================================
RawDecoder::RawDecoder(const vector<uint32_t> &buffer, int start, int end)
{
  mAPVRawSingleEvent.clear();
  mAPVRawHisto.clear();

  fBuf = end - start;
  buf = new unsigned int[fBuf];
  int bufp = 0;
  for(int i=start;i<end;i++)
    {
      buf[bufp] = buffer[i];
      bufp++;
    }
  if(bufp!=fBuf){cout<<"vector passed and vector doesnt match"<<endl;}
  if(fBuf <= 0)
    {
      cout<<"empty vector passed in..."<<endl;
      return;
    }
  Decode();
};

//==========================================================================
RawDecoder::~RawDecoder()
{
  //free buf
  delete[] buf;

  //clear maps
  map<int, vector<int> >::iterator it;
  for(it=mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
  {
      it->second.clear();
  }
  mAPVRawSingleEvent.clear();

  //clear APV raw histos
  map<int, TH1F*>::iterator it_apv;
  for(it_apv=mAPVRawHisto.begin(); it_apv!=mAPVRawHisto.end(); ++it_apv)
  {
      TH1F *h = it_apv->second;
      h->Delete();
  }
  mAPVRawHisto.clear();
  //
}

//==========================================================================
void RawDecoder::Decode()
{
  unsigned int word32bit;

  int mpdid;
  int adc_ch;
  int hybridID;

  //map<int, vector<int> > mpd_event;
  vector<int> apv_event;
  vector<int> apv_margin;
  vector<int> mpd_margin;

  mpd_margin.clear();
  apv_margin.clear();
  
  //find MPD margin, find apv margin

#define BLOCK_HEADER    0x0	// {3'h0, MODULE_ID, EVENT_PER_BLOCK, BLOCK_CNT[7:0]}
#define BLOCK_TRAILER   0x1	// {3'h1, 1'b0, BlockWordCounter}
#define EVENT_HEADER    0x2	// {3'h2, 1'b0, EventCounterFifo_Data}
#define TRIGGER_TIME   0x3	// {3'h3, 1'b0, TimeCounterFifo_Data[39:20]}
  //#define TRIGGER_TIME2   0x3	// {3'h3, 1'b1, TimeCounterFifo_Data[19:0]}
#define APV_CH_DATA     0x4	// {3'h4, ChannelData[20:0]}
#define EVENT_TRAILER   0x5	// {3'h5, 1'b0, LoopDataCounter[11:0], TRIGGER_TIME_FIFO}
#define DATA_NOT_VALID  0x6	// {3'h6, 21'b0}
#define FILLER_WORD     0x7	// {3'h7, 21'b0}

  //find apv margin
    for(int i=0;i<fBuf;i++)//skipping ssp data
  {
    uint32_t data = buf[i];
    uint32_t header;
    uint32_t apv_header;
    header = (data & 0x00e00000)>>21; 
    switch(header)
      {
      case BLOCK_HEADER:
	mpdid=(data&0x001F0000) >> 16;
	//	cout<<"MPDID: "<<mpdid<<endl;
	break;
      case EVENT_HEADER:
	break;
      case TRIGGER_TIME:
	break;
      case APV_CH_DATA:
	switch((data& 0x00180000)>>19)
	  {
	  case 0: //apv header  -------  {1'b0, 1'b0, 1'b0, MEAN[11], DATA_IN[12:0], CH_ID[3:0]};
	    adc_ch=(data&0xf);
	    hybridID=(mpdid<<12)|(adc_ch<<8);
	    //   cout<<"adc_ch: "<<adc_ch<<endl;
	    break;
	  case 1: //data  -------  {1'b0, 1'b1, THRESHOLD_ADDRESS[6:0], data_minus_baseline[11:0]};
	    mAPVRawSingleEvent[hybridID].push_back(data & 0x00000fff);
	    break;
	  case 2: //apv trailer  -------  {1'b1, 1'b0, 2'b0, MODULE_ID[4:0], DATA_IN[11:0]};
	                                                                   //DATA_IN[11:0] = {ApvSampleCounterMinusOne[3:0], frame_counter[7:0]};
	    mAPVRawSingleEvent[hybridID].push_back((data&0xf00)>>8);
	    //	    cout<<" sample count: "<<((data&0xf00)>>8)<<endl;
	    break;
	  case 3: //Trailer  -------  {1'b1, 1'b1, MEAN[10:0], word_count[7:0]};
	    break;
	  default:
	    break;
	  }
	break;
      case EVENT_TRAILER:
	break;
      case BLOCK_TRAILER:
	break;
      case DATA_NOT_VALID:
	break;
      case FILLER_WORD:
	break;
      default:
	break;
      }

  }
}




//==========================================================================
map<int, vector<int> > RawDecoder::GetDecoded()
{
  return mAPVRawSingleEvent;
}

//===========================================================================
map<int, TH1F* > RawDecoder::DrawRawHisto(TCanvas *c)
{
  int mpd_id=0;
  int adc_ch=0;
  int hybridID=0;

  // TCanvas *c = new TCanvas("raw data");
  // c->Divide(12,6);
  map<int, vector<int> >::iterator it;
  for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
    {
      hybridID=it->first;
      mpd_id = GetMPD_ID(hybridID);     
      adc_ch = GetADC_ch(hybridID);
      vector<int> adc_temp = it->second;
      
      int N = adc_temp.size();//cout<<"adc_tempsize:"<<N<<endl;
      //cout<<"mpdid: "<<mpd_id<<"  adcCh: "<<adc_ch<<endl;
      TH1F* h = new TH1F(Form("mpd_%d_ch_%d",mpd_id, adc_ch), Form("mpd_%d_ch_%d_raw_data",mpd_id, adc_ch), 780, 0, 779);
      for(int i=0;i<N;i++)
	{
	  h->Fill(i+1, (Float_t) adc_temp[i]);
	}
      mAPVRawHisto[hybridID] = h;
    }
  map<int, TH1F*>::iterator itRaw;
  int i=1;

  for(itRaw = mAPVRawHisto.begin();itRaw!=mAPVRawHisto.end();itRaw++)
    {
      c->cd(i);
      itRaw->second->SetTitleSize(0.2);
      itRaw->second->SetMaximum(2000);
      itRaw->second->SetMinimum(0);
      itRaw->second->Draw();
      i++;
    }
  c->Update();
 
  int typc = getchar();

  if ((typc==83) || (typc==115)) { // if s or S has pressed, it save the histos in a png file
    c->SaveAs("Result/rawhisto.png");
    getchar(); // return key
  }

  if ((typc==80) || (typc==112)) { // if p or P has pressed, it save the histos in a pdf file
    c->SaveAs("Result/rawhisto.pdf");
    getchar(); // return key
  }

  return mAPVRawHisto;
}

map<int, vector<int> > RawDecoder::GetStripTsAdcMap()
{
  int mpd_id=0;
  int adc_ch=0;
  int hybridID=0;


  map<int, vector<int> >::iterator it;
  for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
  {
    hybridID=it->first;mpd_id= GetMPD_ID(hybridID);
    adc_ch = GetADC_ch(hybridID);
    vector<int> adc_temp = it->second;
    //cout<<adc_temp.size()<<endl;//774
    int TSsize=adc_temp.size()/129;
    //cout<<TSsize<<endl;//774
    for(int i=0; i<TSsize;i++)
      {
	vector<int> singleTSadc_temp;
	singleTSadc_temp.insert(singleTSadc_temp.end(),&adc_temp[129*i],&adc_temp[129*(i+1)]);
	//cout<<"singleTSADCSIZE: "<<singleTSadc_temp.size()<<"     ";
	vector<int> singleTSadc_temp_sorting; 
	//	for(int j=0; j<128;j++)
	//	  {
	//	if(mpd_id==8&&adc_ch==0){cout<<dec<<"i:"<<i<<"singleTSadc_temp:"<<singleTSadc_temp[j]<<endl;}
	//	  }
	singleTSadc_temp_sorting.insert(singleTSadc_temp_sorting.end(),singleTSadc_temp.begin(),singleTSadc_temp.end());
	sort(singleTSadc_temp_sorting.begin(),singleTSadc_temp_sorting.end()-1);
	int iCommonMode=0;
	for ( int k=28; k <100; k++) //omitting largest 4 channels for common mode, necessary to calculate this from the middle
	  {
	    iCommonMode+=singleTSadc_temp_sorting[k];
	    //if((i==0|i==1)&&adc_ch==0)cout <<i<<"  "<<iCommonMode << " ";//if(k==singleTSadc_temp.size()-1){cout<<endl;}
	  }
	iCommonMode = iCommonMode/72; //cout<<"i "<<i<<"commonmode: "<<iCommonMode<<endl;
	for ( int k=0; k <singleTSadc_temp.size()-1; k++) 
	  {
	    singleTSadc_temp[k]-=iCommonMode;
	    //cout<<"commonmode: "<<singleTSadc_temp[k]<<endl;
	  }
	int temphybridID;
	for(int j=0; j<128;j++)
	  {
	    
	    //cout<<"hybridid:"<<hex<<hybridID;
	    temphybridID=hybridID|j;     //if(adc_ch==0){cout<<hex<<"temphybridid:"<<temphybridID<<endl;}
	    mPedestalTsAdc[temphybridID].push_back(singleTSadc_temp[j]);//singleTSadc_temp[j];
	    //cout <<singleTSadc_temp[j]<<endl;
	    //   if(mpd_id==8&&adc_ch==0){cout<<dec<<"i:"<<i<<"singleTSadc_temp:"<<singleTSadc_temp[j]<<endl;}
	    //if((i==0|i==1)&&adc_ch==0)cout <<i<<" XXX  "<<singleTSadc_temp[j] << " ";
	  }
      }
  }// 
  return mPedestalTsAdc;
}

map<int, vector<int> > RawDecoder::ZeroSup(map<int,vector<int> > mMapping, map<int,vector<int> > mPedestalMean, map<int,vector<int> > mPedestalRMS)
{

  map<int, vector<int> > mmHit;

  int mpd_id=0;
  int adc_ch=0;
  int hybridID=0;

  map<int, vector<int> >::iterator it;
  for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
  {
    hybridID = it->first;
    mpd_id = GetMPD_ID(hybridID);     
    adc_ch = GetADC_ch(hybridID);
    
    vector<int> adc_temp = it->second;
    int N = adc_temp.size();//774
    int TSsize=N/129;
    //cout<<TSsize<<endl;//774
    int CommonMode[TSsize];


    for(int i=0;i<N;i++)
	{
	  if((i%129)!=128) adc_temp[i]-=mPedestalMean[hybridID][i%129];
	}


    for(int i=0; i<TSsize;i++)
      { 
	CommonMode[i]=0;
	vector<int> singleTSadc_temp_sorting;
	singleTSadc_temp_sorting.insert(singleTSadc_temp_sorting.end(),&adc_temp[129*i],&adc_temp[129*(i+1)]);
	  	 
	//cout<<(singleTSadc_temp_sorting.size()-5)<<endl;
	sort(singleTSadc_temp_sorting.begin(),singleTSadc_temp_sorting.end()-1);

	for ( int k=28; k <100; k++) //omitting largest 4 channels for common mode, necessary to calculate this from the middle
	  {  
	    CommonMode[i]+=singleTSadc_temp_sorting[k];
	    //CommonMode[i]+=500;

	  }

	CommonMode[i] = CommonMode[i]/72; 
	//	cout<<"Commonmode: "<<CommonMode[i];
      }
 
    for(int j=0; j<128;j++)
      {
	int adcSum_temp=0;
	for(int i=0;i<TSsize;i++)
	  {//cout<<"ADC:"<<adc_temp[j+129*i]<<"    ";
	    adcSum_temp = adcSum_temp+adc_temp[j+129*i]-CommonMode[i];
	    //cout<<"ADC:"<<adcSum_temp<<"    ";
	  }
	// cout<<endl<<endl;
	adcSum_temp = adcSum_temp/TSsize; 
	 
	//cout<<j<<"Mean  "<<mPedestalMean[mpd_id][adc_ch][j]<<"  "<<mPedestalRMS[mpd_id][adc_ch][j]<<endl;

	if(adcSum_temp>5*mPedestalRMS[hybridID][j])
	  { 
	  
	    int RstripPos=j;	  
	    int RstripNb = ChNb[j];
	    //int RstripNb=32*(j%4)+8*(int)(j/4)-31*(int)(j/16);                        //channel re-matching for apv25 chip
	    ////stripNb=(8*(int)(stripNb/4)+3-stripNb)*((int)(stripNb/4)%2)+stripNb*(1-((int)(stripNb/4)%2));
	    //RstripNb=RstripNb+1+RstripNb%4-5*(((int)(RstripNb/4))%2);                 //channel re-matching for INFN type APV front-end card
	    //cout<<RstripNb<<", ";
	    RstripNb=RstripNb+(127-2*RstripNb)*mMapping[hybridID][3];                   //re-matching for inverted strips Nb
	    RstripPos=RstripNb+128*mMapping[hybridID][2];                               // calculate position

	    int detID = mMapping[hybridID][0];
	    int planeID = mMapping[hybridID][1];
	    int HitHybridID = (detID<<13)|(planeID<<12)|RstripPos;
	    for(int i=0; i<TSsize;i++)
	      { 
		mmHit[HitHybridID].push_back(adc_temp[j+129*i]-CommonMode[i]);//if stop here, reduce 200ms/10k event
		//	cout<<(adc_temp[j+129*i]-CommonMode[i]-mPedestalMean[hybridID][j])<<"  ";
	      }     
	    //	    cout<<endl;
	    //  j++;
	  } 

      }  
  }// 
  return mmHit;
}

map<int, vector<int> > RawDecoder::DrawHits(map<int,vector<int> > mMapping, map<int,vector<int> > mPedestalMean, map<int,vector<int> > mPedestalRMS, TCanvas *c)
{
  map<int, vector<int> > mmHit;

  int mpd_id=0;
  int adc_ch=0;
  int hybridID=0;

  map<int,TH1F*> hitHisto;
  for(int i=0;i<10;i++)
    {
      //    hitHisto[i]->Delete();
      hitHisto[i]= new TH1F(Form("GEM_%d__plane_%d",i/2,i%2),Form("GEM_%d__plane_%d",i/2,i%2),1600,0,1600);
    }



  map<int, vector<int> >::iterator it;
  for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
  {
    hybridID = it->first;
    mpd_id = GetMPD_ID(hybridID);     
    adc_ch = GetADC_ch(hybridID);
    
    vector<int> adc_temp = it->second;
    int N = adc_temp.size();//774
    int TSsize=N/129;
    //cout<<TSsize<<endl;//774
    int CommonMode[TSsize];


    for(int i=0;i<N;i++)
	{
	  if((i%129)!=128) adc_temp[i]-=mPedestalMean[hybridID][i%129];
	}


    for(int i=0; i<TSsize;i++)
      { 
	CommonMode[i]=0;
	vector<int> singleTSadc_temp_sorting;
	singleTSadc_temp_sorting.insert(singleTSadc_temp_sorting.end(),&adc_temp[129*i],&adc_temp[129*(i+1)]);
	  	 
	//cout<<(singleTSadc_temp_sorting.size()-5)<<endl;
	sort(singleTSadc_temp_sorting.begin(),singleTSadc_temp_sorting.end()-1);

	for ( int k=10; k <118; k++) //omitting largest 4 channels for common mode, necessary to calculate this from the middle
	  {  
	    CommonMode[i]+=singleTSadc_temp_sorting[k];
	    //CommonMode[i]+=500;

	  }

	CommonMode[i] = CommonMode[i]/108; 
	//	cout<<"Commonmode: "<<CommonMode[i];
      }
 
    for(int j=0; j<128;j++)
      {
	int adcSum_temp=0;
	for(int i=0;i<TSsize;i++)
	  {//cout<<"ADC:"<<adc_temp[j+129*i]<<"    ";
	    adcSum_temp = adcSum_temp+adc_temp[j+129*i]-CommonMode[i];
	    //cout<<"ADC:"<<adcSum_temp<<"    ";
	  }
	// cout<<endl<<endl;
	adcSum_temp = adcSum_temp/TSsize; 
	 
	//cout<<j<<"Mean  "<<mPedestalMean[mpd_id][adc_ch][j]<<"  "<<mPedestalRMS[mpd_id][adc_ch][j]<<endl;

	if(adcSum_temp>6*mPedestalRMS[hybridID][j])
	  { 
	  
	    int RstripPos=j;	  
	    int RstripNb = ChNb[j];
	    //int RstripNb=32*(j%4)+8*(int)(j/4)-31*(int)(j/16);                        //channel re-matching for apv25 chip
	    ////stripNb=(8*(int)(stripNb/4)+3-stripNb)*((int)(stripNb/4)%2)+stripNb*(1-((int)(stripNb/4)%2));
	    //RstripNb=RstripNb+1+RstripNb%4-5*(((int)(RstripNb/4))%2);                 //channel re-matching for INFN type APV front-end card
	    //cout<<RstripNb<<", ";
	    RstripNb=RstripNb+(127-2*RstripNb)*mMapping[hybridID][3];                   //re-matching for inverted strips Nb
	    RstripPos=RstripNb+128*mMapping[hybridID][2];                               // calculate position

	    int detID = mMapping[hybridID][0];
	    int planeID = mMapping[hybridID][1];
	    int HitHybridID = (detID<<13)|(planeID<<12)|RstripPos;

	   
	    hitHisto[detID*2+planeID]->Fill(RstripPos,adcSum_temp);
	    //  j++;
	  } 

      }  
  }// 

  for(int i=0;i<6;i++)
    {
      c->cd(i+1);
      hitHisto[i]->Draw();
    }
  c->Update();
  getchar();
    for(int i=0;i<10;i++)
    {
      hitHisto[i]->Delete();
    }

  return mmHit;
}
