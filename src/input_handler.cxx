#include "input_handler.h"
#include "GEMConfigure.h"
#include <stdint.h>
#include <fstream>
#include <TStyle.h>

#include <cmath>

//_______________________________________________
InputHandler::InputHandler()
{
  Initialize();
  fRawDecoder = 0;
}

//______________________________________________
InputHandler::~InputHandler()
{
/*
  delete[] Vstrip;
  delete[] VdetID;
  delete[] VplaneID;
  delete[] adc0 ;
  delete[] adc1 ;
  delete[] adc2 ;
  delete[] adc3 ;
  delete[] adc4 ;
  delete[] adc5 ;
*/
}

//##########################################################################################################################
//######################################################  Initialization  ##################################################
//##########################################################################################################################

//______________________________________________
int InputHandler::Initialize()
{
  cout<<"Initializing!!!"<<endl;
  LoadConfigFile();
  LoadMapping();

  cout<<"runtype = "<< runtype << endl;
  if(runtype=="CalPedestal")
    {
      InitPedestalToWrite();
      return 1;
    }
  if(runtype=="GetRootFile") 
    {
      LoadPedestalToRead();
      InitRootFile();
      return 1;
    }
  if(runtype=="RawDataMonitor") 
    {
      gStyle->SetOptStat(0);
      cRaw = new TCanvas("Raw_frame_Monitor","Raw Frame Monitor",0,0,2000, 2000);
      cRaw->Divide(raw_co,raw_row);
      //cRaw->Divide(15,16,0.000001, 0.000001);
      /*  
      cRaw0 = new TCanvas("MPD_0","MPD 0",0,0,1600,800);
      cRaw0->Divide(4,4);

      cRaw1 = new TCanvas("MPD_1","MPD 1",0,0,1600,800);
      cRaw1->Divide(4,4);

      cRaw2 = new TCanvas("MPD_2","MPD 2",0,0,1600,800); 
      cRaw2->Divide(4,4);
     
      cRaw3 = new TCanvas("MPD_3","MPD 3",0,0,1600,800);
      cRaw3->Divide(4,4);
      
      cRaw4 = new TCanvas("MPD_4","MPD 4",0,0,1600,800);
      cRaw4->Divide(4,4);

      cRaw5 = new TCanvas("MPD_5","MPD 5",0,0,1600,800);
      cRaw5->Divide(4,4);

      cRaw6 = new TCanvas("MPD_6","MPD 6",0,0,1600,800);
      cRaw6->Divide(4,4);

      cRaw7 = new TCanvas("MPD_7","MPD 7",0,0,1600,800);
      cRaw7->Divide(4,4);

      cRaw8 = new TCanvas("MPD_8","MPD 8",0,0,1600,800);
      cRaw8->Divide(4,4);

      cRaw9 = new TCanvas("MPD_9","MPD 9",0,0,1600,800);
      cRaw9->Divide(4,4);

      cRaw10 = new TCanvas("MPD_10","MPD 10",0,0,1600,800);
      cRaw10->Divide(4,4);

      cRaw11 = new TCanvas("MPD_11","MPD 11",0,0,1600,800);
      cRaw11->Divide(4,4);
      */
    }
  if(runtype=="singleEventHitMonitor") 
    {
      LoadPedestalToRead();
      cHit = new TCanvas("Hit_Monitor","Hit Monitor",0,0,1000,800);
      cHit->Divide(2,3);
    }
  return 1;
}


//______________________________________________
int InputHandler::LoadConfigFile()
{
  GEMConfigure *configure = new GEMConfigure;
  configure->LoadConfigure();
  runtype = configure->GetRunType();
  if(runtype=="RawDataMonitor"){raw_co = configure->GetRawCo();raw_row = configure->GetRawRow();}
  PedFiletoSAVE = configure->GetSavePedPath();
  PedFiletoLOAD = configure->GetLoadPedPath();
  HitRootFiletoSAVE=configure->hitrootfile_path;
  MappingFile = configure->GetMapping();
  nEvt=configure->GetNumEvt();
  NbFile=configure->nFile;
  for(int i=0;i<NbFile;i++)
    {  
      filenameList[i]=configure->fileList[i];
      NbofInputFile++;
    }
  cout<<"About to decode "<<NbFile<<" files:  "<<endl;
  cout<<configure->fileHeader<<"."<<configure->evioStart<<" to "
      <<configure->fileHeader<<"."<<configure->evioEnd<<endl;
  return 1;
}

//______________________________________________
int InputHandler::LoadMapping()
{
  
  char *MappingFile_temp = new char[100];
  strcpy(MappingFile_temp,MappingFile.c_str());
  cout<<"Mapping!!!   "<<MappingFile_temp<<endl;
  ifstream filestream (MappingFile_temp, ifstream::in); 
  cout<<filestream.is_open()<<endl;
  delete[] MappingFile_temp;
  string line;
  int Mapping_mpdId,
      Mapping_ADCId,
      Mapping_I2C,
      Mapping_GEMId,
      Mapping_Xis,
      Mapping_Pos,
      Mapping_Invert,
      Mapping_hybridID;
  while (getline(filestream, line) ) {
    cout << "Ingest Map Line: " << line << endl;
    line.erase( remove_if(line.begin(), line.end(), ::isspace), line.end() );
    if( line.find("#") == 0 ) continue; 
    char *tokens = strtok( (char *)line.data(), ",");
    if(tokens !=NULL){
      cout<<tokens<<" ";Mapping_mpdId=atoi(tokens);
      tokens = strtok(NULL, " ,");cout<<tokens<<" ";Mapping_GEMId=atoi(tokens);
      tokens = strtok(NULL, " ,");cout<<tokens<<" ";Mapping_Xis=atoi(tokens);
      tokens = strtok(NULL, " ,");cout<<tokens<<" ";Mapping_ADCId=atoi(tokens);
      tokens = strtok(NULL, " ,");cout<<tokens<<" ";Mapping_I2C=atoi(tokens);
      tokens = strtok(NULL, " ,");cout<<tokens<<" ";Mapping_Pos=atoi(tokens);
      tokens = strtok(NULL, " ,");cout<<tokens<<" ";Mapping_Invert=atoi(tokens);
      Mapping_hybridID = (Mapping_mpdId<<12)|(Mapping_ADCId<<8);
      mMapping[Mapping_hybridID].push_back(Mapping_GEMId);//0
      mMapping[Mapping_hybridID].push_back(Mapping_Xis);//1
      mMapping[Mapping_hybridID].push_back(Mapping_Pos);//2
      mMapping[Mapping_hybridID].push_back(Mapping_Invert);//3
      //cout<<"test: "<<mMapping[Mapping_mpdId][Mapping_ADCId][2]<<endl;
    }
    cout<<endl;
  }
  filestream.close();
  cout << "Map loaded" << endl;
  return 1;
}

//______________________________________________
int InputHandler::InitPedestalToWrite()
{
   for(map<int,vector<int> >::iterator it = mMapping.begin(); it!=mMapping.end(); ++it)
	{
	  int hybridID = it->first;
	  int mpd_id = GetMPD_ID(hybridID);
	  int adc_ch = GetADC_ch(hybridID);
	  for(int stripNb=0;stripNb<128;stripNb++)
	    {
	      mPedestalHisto[hybridID|stripNb] = new TH1F(Form("mpd_%d_ch_%d_Strip_%d",mpd_id, adc_ch,stripNb), Form("mpd_%d_ch_%d_Strip_%d_pedestal_data",mpd_id, adc_ch,stripNb), 3500, -500, 3000);
	    }
	  mPedestalMean[hybridID] = new TH1F(Form("PedestalMean(offset)_mpd_%d_ch_%d",mpd_id, adc_ch), Form("PedestalMean(offset)_mpd_%d_ch_%d",mpd_id, adc_ch), 128, 0, 128);
	  mPedestalRMS[hybridID] = new TH1F(Form("PedestalRMS_mpd_%d_ch_%d",mpd_id, adc_ch), Form("PedestalRMS_mpd_%d_ch_%d",mpd_id, adc_ch), 128, 0, 128);
	      mPedestalRMS[hybridID]->SetMinimum(0);
	}
}

//______________________________________________
int InputHandler::LoadPedestalToRead()
{
  char *PedFilename_temp = new char[100];
  strcpy(PedFilename_temp,PedFiletoLOAD.c_str());
  f = new TFile(PedFilename_temp,"READ" );
  delete[] PedFilename_temp;

  int hybridID,mpd_id,adc_ch;  
  for(map<int,vector<int> >::iterator it = mMapping.begin(); it!=mMapping.end(); ++it)
	{
	  hybridID = it->first;
	  mpd_id = GetMPD_ID(hybridID);
	  adc_ch = GetADC_ch(hybridID);
	  cout<<"loading pedestal from MPD id:"<<mpd_id<<"ADC Ch:"<<adc_ch<<endl;
	  //loading pedestal information from root files
	  TH1F* hMean = (TH1F*)f->Get(Form("PedestalMean(offset)_mpd_%d_ch_%d",mpd_id, adc_ch));
	  TH1F* hRMS  = (TH1F*)f->Get(Form("PedestalRMS_mpd_%d_ch_%d",mpd_id, adc_ch));
	  for(int i=0;i<128;i++)
	    {
	      hMean->GetBinContent(i+1);
	      hRMS->GetBinContent(i+1);
	      mvPedestalMean[hybridID].push_back(hMean->GetBinContent(i+1));
	      mvPedestalRMS[hybridID].push_back(hRMS->GetBinContent(i+1));
	      //  cout<<i<<"Mean  "<<hMean->GetBinContent(i+1)<<"  "<<hRMS->GetBinContent(i+1)<<endl;
	    }
	}
  cout<<"finished loading pedestal"<<endl;
  f->Close();
  //end of loading
}

//______________________________________________
int InputHandler::InitRootFile()
{
  Vstrip   = new Int_t[20000];
  VdetID   = new Int_t[20000];
  VplaneID = new Int_t[20000];
  adc0     = new Int_t[20000];
  adc1     = new Int_t[20000];
  adc2     = new Int_t[20000];
  adc3     = new Int_t[20000];
  adc4     = new Int_t[20000];
  adc5     = new Int_t[20000];
  fCH      = new Int_t[100];
  ftiming  = new Int_t[100];
  fadc     = new Int_t[100];
  tCH      = new Int_t[100];
  ttiming  = new Double_t [100];

  Hit      = new TTree("GEMHit","Hit list");

  //-------------GEM----------------------
  Hit->Branch("evtID",&EvtID,"evtID/I");
  Hit->Branch("nch",&nch,"nch/I");
  Hit->Branch("strip",Vstrip,"strip[nch]/I");
  Hit->Branch("detID",VdetID,"detID[nch]/I");
  Hit->Branch("planeID",VplaneID,"planeID[nch]/I");
  Hit->Branch("adc0",adc0,"adc0[nch]/I");
  Hit->Branch("adc1",adc1,"adc1[nch]/I");
  Hit->Branch("adc2",adc2,"adc2[nch]/I");
  Hit->Branch("adc3",adc3,"adc3[nch]/I");
  Hit->Branch("adc4",adc4,"adc4[nch]/I");
  Hit->Branch("adc5",adc5,"adc5[nch]/I");
  //--------------fadc---------------------
  Hit->Branch("nfadc",&nfadc,"nfadc/I");//nb of fadc channel fired
  Hit->Branch("fCH",fCH,"fCH[nfadc]/I");//slot+ch
  Hit->Branch("ftimting",ftiming,"ftiming[nfadc]/I");//timing
  Hit->Branch("fadc",fadc,"fadc[nfadc]/I");//adc total/adc max
  //--------------tdc----------------------
  Hit->Branch("ntdc",&ntdc,"ntdc/I");//nb of tdc channel fired
  Hit->Branch("tCH",tCH,"tCH[ntdc]/I");//channel nb
  Hit->Branch("ttiming",ttiming,"ttiming[ntdc]/D");//tdc value
}


//##########################################################################################################################
//##################################################  Event Processing  ####################################################
//##########################################################################################################################
//______________________________________________
int InputHandler::ProcessAllFiles()
{ 
  const int mpd_tag=3561,fadc_tag=3,tdc_tag=6;
  int entry = 0;
  for(int NFile=0;NFile<NbofInputFile;NFile++)
    {
      string filename = filenameList[NFile];
      cout<<"Processing "<<filename<<endl;
      try{
	evioFileChannel chan(filename.c_str(), "r");
	chan.open();
	cout<<">"<<endl;

	while(chan.read()&&entry<nEvt)
	  {
	    
	    map<int, Double_t > mTDC;
	    map<int, vector<int> > mFADc;//slot+channel
	    evioDOMTree event(chan);
	    //  cout<<event.toString()<<endl;
	    // cout<<"##"<<endl;
	    // getchar();
	    evioDOMNodeListP mpdEventList = event.getNodeList( isLeaf() );
	    //evioDOMNodeListP pList = event.getNodeList( );
	    //for_each(pList->begin(), pList->end(), toCout());
	    //cout<<"##"<<endl;
	    //getchar();

	    if(EvtID%100>=0)
	      {
		cout<<"event: "<<EvtID<<endl;
		cout<<"number of banks: "<<mpdEventList->size()<<endl;
	      }
	    evioDOMNodeList::iterator iter;
	    for(iter=mpdEventList->begin(); iter!=mpdEventList->end(); ++iter)
	      {
		int banktag = (*iter)->tag;
		
       		switch(banktag){
		case mpd_tag:
		  vector<uint32_t> *block_vec = (*iter)->getVector<uint32_t>();
		  entry+=ProcessSspBlock(*block_vec);
		  break;
		  //	default: 
		  // break;
		}
	      }
	  
	    if(runtype=="GetRootFile") 
	      {	
		Hit->Fill();
	      }
	  }
	chan.close();

      } catch (evioException e) 
	{
	  cerr <<endl <<e.toString() <<endl <<endl;
	  exit(EXIT_FAILURE);
	}
    }
  //exit(EXIT_SUCCESS);
  return entry;
}

//______________________________________________
int InputHandler::ProcessSspBlock(const vector<uint32_t> &block_vec)
{
  int vec_size = block_vec.size();
  if(EvtID%100==0) cout<<"procBlock: blockvec size: "<<vec_size<<endl;
  int NofEvtInBlk=0;

  NofEvtInBlk++;
  ProcessSingleSspEvent(block_vec,0);
  EvtID++;

  // ssp multiple event block
  /*for(int istart=0; istart<vec_size; istart++)
    {
      if(((block_vec[istart]>>24)&0xf8)==0x90)// event header
	{
	  NofEvtInBlk++;
	  istart = ProcessSingleSspEvent(block_vec,istart);
	  EvtID++;
	}
		    
    }
  */
  return NofEvtInBlk;
}

//______________________________________________
int InputHandler::ProcessSingleSspEvent(const vector<uint32_t> &block_vec, int istart)
{
  int iend;
  int vec_size = block_vec.size();
  for(iend=istart+1;iend<vec_size;iend++)
    {
      uint32_t tag = (block_vec[iend]>>24)&0xf8;
      if(tag==0x90||tag==0x88) break;
    }  
  RawDecoder raw_decoder(block_vec,istart,iend);
  
  //if(runtype=="RawDataMonitor")  {raw_decoder.DrawRawHisto(cRaw, cRaw0,cRaw1,cRaw2,cRaw3,cRaw4,cRaw5,cRaw6,cRaw7,cRaw8,cRaw9,cRaw10,cRaw11);}
  if(runtype=="RawDataMonitor")  {raw_decoder.DrawRawHisto(cRaw);}
  if(runtype=="CalPedestal")     
    {
      map<int, vector<int> > mTsAdc = raw_decoder.GetStripTsAdcMap();
      CalSinglePedestal(mTsAdc);
    }
  if(runtype=="GetRootFile")   
    {
      map<int, vector<int> > mmHit = raw_decoder.ZeroSup(mMapping,mvPedestalMean,mvPedestalRMS);
      FillRootTree(mmHit);
    }
  if(runtype=="singleEventHitMonitor")
    {
      raw_decoder.DrawHits(mMapping,mvPedestalMean,mvPedestalRMS,cHit);
      
      // cHit->Update();
    }


  //entry++;

  return (iend-1);
}

//______________________________________________
int InputHandler::CalSinglePedestal(const map<int, vector<int> > &mTsAdc)
{
  //Fill strip average ADC Histogram for calculating pedestal
  for( map<int,vector<int> >::const_iterator it = mTsAdc.begin(); it!=mTsAdc.end(); ++it)
    {
      int hybridID = it->first;
      vector<int> adc_temp = it->second;
      int adcSum_temp=0;
      int TSsize = adc_temp.size(); 	  
      for(int i=0; i<TSsize;i++)
	{ 
	  adcSum_temp+=adc_temp[i];
	}
      adcSum_temp=adcSum_temp/TSsize;
      mPedestalHisto[hybridID]->Fill(adcSum_temp);
    }
  return 1;
}

//______________________________________________
int InputHandler::FillRootTree(const map<int, vector<int> > &mmHit)
{
  vector<int> stripVector;
  Int_t nstrip=0;
   
  for(map<int, vector<int> >::const_iterator it = mmHit.begin(); it!=mmHit.end(); ++it)
    {
      //cout<<"planeid: "<<planeid<<"  "<<stripVector.size()<<"DDDDD"<<endl;
      //cout<<"HHH"<<ittt->first<<endl;
      stripVector.push_back(it->first);
      // cout<<EvtID<<endl;
      //getchar();
    }

  for(int i=0;i<stripVector.size();i++)
    {
      Vstrip[nstrip]=GetStripPos(stripVector[i]);
      adc0[nstrip]=mmHit.at(stripVector[i])[0];
      adc1[nstrip]=mmHit.at(stripVector[i])[1];
      adc2[nstrip]=mmHit.at(stripVector[i])[2];

      adc3[nstrip]=mmHit.at(stripVector[i])[3];
      adc4[nstrip]=mmHit.at(stripVector[i])[4];
      adc5[nstrip]=mmHit.at(stripVector[i])[5];
      //cout<<"ADC:"<<stripVector[i]<<" "<<adc0[nstrip]<<" "<<adc1[nstrip]<<" "<<adc2[nstrip]<<" "<<adc3[nstrip]<<" "<<adc4[nstrip]<<" "<<adc5[nstrip]<<endl;
      VdetID[nstrip]=GetDet_ID(stripVector[i]); //cout<<VdetID[nstrip]<<endl;
      VplaneID[nstrip]=GetPlane_ID(stripVector[i]); //cout<<"planeid"<<VplaneID[nstrip]<<endl;
	 
      //	  for(int j=0;j<6;j++)
      // {
      //  txt<<setw(12)<<setfill(' ')<<entry<<setw(12)<<setfill(' ')<<(VdetID[nstrip]*2+VplaneID[nstrip])<<setw(12)<<setfill(' ')<<Vstrip[nstrip]<<setw(12)<<setfill(' ')<<mmHit[stripVector[i]][j]<<endl;
      //}


      //	if(CountFlag[(2*VdetID[nstrip]+VplaneID[nstrip])]!=0){hitcount[(2*VdetID[nstrip]+VplaneID[nstrip])]+=1;}
      //	CountFlag[(2*VdetID[nstrip]+VplaneID[nstrip])]=0;
	  
      nstrip++;
    }
  nch=nstrip;
  //cout<<"nch:: "<<nch<<endl;
  //if(nch==0) return 1;
  //EvtID=entry;
  //cout<<EvtID<<" :: "<<endl;
  //Hit->Fill();
  return 1;
}



//##########################################################################################################################
//#####################################################  Post Processing  ##################################################
//##########################################################################################################################
//______________________________________________
int InputHandler::Summary()
{
  if(runtype=="CalPedestal") SummPedestal();
  if(runtype=="GetRootFile") SummRootFile();
  return 1;
}

//______________________________________________
int InputHandler::SummPedestal()
{

  //Calculating pedestal
  char *PedFilename_temp = new char[100];
  std::strcpy(PedFilename_temp,PedFiletoSAVE.c_str());
  f = new TFile(PedFilename_temp,"RECREATE");
  delete[] PedFilename_temp;
  map<int,TH1F* > ::iterator it;
  for(it = mPedestalHisto.begin(); it!=mPedestalHisto.end(); ++it)
    {
      int hybridID = it->first;//cout<<"hybridid:"<<hybridID<<endl;
      int mpd_id = GetMPD_ID(hybridID);// cout<<"mpdid:"<<mpd_id<<endl;
      int adc_ch = GetADC_ch(hybridID);// cout<<"adcid:"<<adc_ch<<endl;
      int stripNb = GetStripNb(hybridID); //if(adc_ch==0){cout<<"stripid:"<<stripNb<<endl;}     
      TH1F* Pedestal_temp = it->second;
      float mean = Pedestal_temp->GetMean();//if(adc_ch==0){cout<<"Mean:  "<<mean<<endl;}
      float rms  = Pedestal_temp->GetRMS();
      //cout<<"StripNb: "<<stripNb<<" RMS: "<<rms<<"MEAN: "<<mean<<endl;
      mPedestalMean[hybridID&0xfff00]->Fill(stripNb,mean);
      mPedestalRMS[hybridID&0xfff00]->Fill(stripNb,rms);      
    }
  
   for(map<int,vector<int> >::iterator it = mMapping.begin(); it!=mMapping.end(); ++it)
	{
	  int hybridID = it->first;
	  mPedestalMean[hybridID]->SetDirectory(f);
	  mPedestalRMS[hybridID]->SetDirectory(f);
	}

  f->Write();
  f->Close();
  //end of Calculating pedestal

  //delete histograms
  for(map<int,vector<int> > ::iterator it = mTsAdc.begin(); it!=mTsAdc.end(); ++it)
    {
      mPedestalHisto[it->first]->Delete();
    }
}

//______________________________________________
int InputHandler::SummRootFile()
{ 
  char *HitFilename_temp = new char[100];
  std::strcpy(HitFilename_temp,HitRootFiletoSAVE.c_str());
  TFile *Hit_rootfile = new TFile(HitFilename_temp,"RECREATE");
  delete[] HitFilename_temp;

  Hit->SetDirectory(Hit_rootfile);
  cout<<"writing root file..."<<endl;
  Hit->Write();
  Hit_rootfile->Write();
  Hit_rootfile->Close();
}
