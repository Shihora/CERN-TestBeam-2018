//root
#include <TLine.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <THStack.h>
#include <THistPainter.h>
//#include <TStyle.h>

//C, C++
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
//#include <stdlib.h>
//#include <string>
//#include <iomanip>

//specific
#include "geometry.h"
#include "analysis.h"
#include "read.h"

float SP = 0.3125;
float pe = 47.46;//mV*ns
//vector<float> pe_SiPM = {32.14, 39.33, 34.20, 30.79, 34.09, 29.99, 30.69, 29.95}; //a,b,c,d,e,f,g,h  -  Gain-Baseline from fit
vector<float> pe_SiPM = {42.01, 34.67, 34.28, 33.84, 37.55, 34.68, 33.81, 38.84}; //sorted by Wavecatcher-Channel
vector<float> SiPM_shift = {2.679, 2.532, 3.594, 3.855, 3.354, 3.886, 3.865, 4.754};
int wavesPrintRate = 100;
int ch0PrintRate = 1000000;
int trigPrintRate = 1000000;//100
int signalPrintRate = 100000;//100
double coef = 2.5 / (4096 * 10);
string WCHU ("WC-Version:1.14"), WCAlexander ("WC-Version:1.7");

extern string WCVersion;
extern int runNr;
extern float horizontal;
extern float vertical;
extern float angle;
extern int pdgID;
extern float energy;
extern int isSP;
extern int mp;
extern int safPMT2;
extern int safPMT1;
extern int safSiPM;
extern int trackL;

void read(TString _inFileList, TString _inDataFolder, TString _outFile){

  TF1* fTrigFit = new TF1("fTrigFit","gaus");
  fTrigFit->SetParameter(0,800);
  fTrigFit->SetParameter(2,1);
  fTrigFit->SetLineWidth(1);
  
  
  
  ///////////////////Root file with data/////////////////
  TFile *rootFile = new TFile( _outFile, "RECREATE");
  if (rootFile->IsZombie()) {
    cout << "PROBLEM with the initialization of the output ROOT ntuple "
	 << _outFile << ": check that the path is correct!!!"
	 << endl;
    exit(-1);
  }
  TTree *tree = new TTree("T", "USBWC Data Tree");
  //rootFile->SetCompressionLevel(2);
  //tree->SetAutoSave(1000000);
  // Create new event
  TTree::SetBranchStyle(0);

  Int_t EventNumber=-999;
  Int_t LastEventNumber=-999;
  Float_t SamplingPeriod = -999;
  Double_t EpochTime = -999;
  Int_t Year = -999;
  Int_t Month = -999;
  Int_t Day = -999;
  Int_t Hour = -999;
  Int_t Minute = -999;
  Int_t Second = -999;
  Int_t Millisecond = -999;
  Float_t trigT = -999;//t_trig = (t0+t1+t2+t3)/4
  Float_t tPMT1 = -999;
  Float_t tPMT2 = -999;
  Float_t tPMT2i = -999;
  Float_t tSUMp = -999;
  Float_t tSUMm = -999;
  Float_t trigTp = -999;//t_trig' = [(t0+t1)-(t2+t3)]/4
  Float_t t0t1 = -999;//t0t1 = [(t0-t1)]
  Float_t t2t3 = -999;//t2t3 = [(t2-t3)]
  Int_t isVeto = -999; //variable to define veto, 1 if veto, 0 if not, -999 if undefined
  Int_t isTrig = -999;
  Int_t isLastEvt = -999;
  Int_t isGoodSignal_5 = -999;
  Float_t trigGate = -999;
  Int_t nCh = -1;
  int nActiveCh = -1;
  Int_t ChannelNr[16];
  std::vector<float> amp(16,-999);
  std::vector<float> max(16,-999);
  std::vector<float> min(16,-999);
  Float_t t[16];
  Float_t tSiPM[16];
  Float_t BL[16];//store baseline for 16 channels
  Float_t BL_RMS[16];//store rms of baseline for 16 channels
  float BL_output[2];//array used for output getBL-function
  float Integral_0_300[16];//array used to store Integral of signal from 0 to 300ns
  float Integral[16];
  float Integral_Correction;
  float Integral_mVns[16];
  int NumberOfBins;
  Int_t EventIDsamIndex[16];
  Int_t FirstCellToPlotsamIndex[16];
  
  std::vector<TH1F*> hChSum;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChSum_%d",i);
    TH1F* h = new TH1F("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h->SetName(name);
    hChSum.push_back(h);
  }
  std::vector<TH1F*> hChShift;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChShift_%d",i);
    TH1F* h = new TH1F("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h->SetName(name);
    hChShift.push_back(h);
  }
  std::vector<TH1F> hChtemp;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChtemp_%d",i);
    TH1F h("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h.SetName(name);
    hChtemp.push_back(h);
  }
  std::vector<TH1F> hChShift_temp;
  for(int i=0;i<16;i++){
    TString name("");
    name.Form("hChShift_temp_%d",i);
    TH1F h("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h.SetName(name);
    hChShift_temp.push_back(h);
  }
  
  Short_t amplValues[16][1024];
  TH1F hCh("hCh","dummy;ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
  TString plotSaveFolder  = _outFile;
  plotSaveFolder.ReplaceAll("out.root","");
  TCanvas cWaves("cWaves","cWaves",1000,700);
  cWaves.Divide(4,4);
  TCanvas cCh0("cCh0","cCh0",1500,900);
  cCh0.Divide(2,2);
  TCanvas cTrig("cTrig","cTrig",1500,900);
  cTrig.Divide(2,2);
  TCanvas cSignal("cSignal","cSignal",1500,900);
  cSignal.Divide(2,2);

  tree->Branch("EventNumber",&EventNumber, "EventNumber/I");
  tree->Branch("SamplingPeriod", &SamplingPeriod,  "SamplingPeriod/F");
  tree->Branch("EpochTime",&EpochTime, "EpochTime/D");
  tree->Branch("Year",&Year, "Year/I");
  tree->Branch("Month",&Month, "Month/I");
  tree->Branch("Day",&Day, "Day/I");
  tree->Branch("Hour",&Hour, "Hour/I");
  tree->Branch("Minute",&Minute, "Minute/I");
  tree->Branch("Second",&Second, "Second_/I");
  tree->Branch("Millisecond",&Millisecond, "Millisecond/I");
  tree->Branch("trigT",&trigT, "trigT/F");
  tree->Branch("tPMT1",&tPMT1, "tPMT1/F");
  tree->Branch("tPMT2",&tPMT2, "tPMT2/F");
  tree->Branch("tPMT2i",&tPMT2i, "tPMT2i/F");
  tree->Branch("tSUMp",&tSUMp, "tSUMp/F");
  tree->Branch("tSUMm",&tSUMm, "tSUMm/F");

  tree->Branch("runNr",&runNr, "runNr/I");//run number in google table
  tree->Branch("horiz",&horizontal,"horiz/F");// horizontal position of the box units: [cm]
  tree->Branch("vert",&vertical,"vert/F");//vertical position of the box, units: [cm]
  tree->Branch("angle",&angle,"angle/F");
  tree->Branch("pdgID",&pdgID,"pdgID/I");
  tree->Branch("energy",&energy,"energy/F");
  tree->Branch("isSP",&isSP,"isSP/I");
  tree->Branch("mp",&mp,"mp/I");
  tree->Branch("safPMT2",&safPMT2,"safPMT2/I");//solid angle factor
  tree->Branch("safPMT1",&safPMT1,"safPMT1/I");//solid angle factor
  tree->Branch("safSiPM",&safSiPM,"safSiPM/I");//solid angle factor
  tree->Branch("trackL",&trackL,"trackL/I");//track length
  tree->Branch("isLastEvt",&isLastEvt,"isLastEvt/I");
  
  tree->Branch("trigGate",&trigGate,"trigGate/F");
  tree->Branch("trigTp",&trigTp, "trigTp/F");
  tree->Branch("t0t1",&t0t1, "t0t1/F");//t0t1 = [(t0-t1)]
  tree->Branch("t2t3",&t2t3, "t2t3/F");
  tree->Branch("isVeto",&isVeto,"isVeto/I");
  tree->Branch("isTrig",&isTrig,"isTrig/I");
  tree->Branch("isGoodSignal_5",&isGoodSignal_5,"isGoodSignal_5/I");
  
  tree->Branch("nCh",&nCh, "nCh/I");
  tree->Branch("ch",ChannelNr, "ch[nCh]/I");
  tree->Branch("amp",amp.data(), "amp[nCh]/F");
  tree->Branch("max",max.data(), "max[nCh]/F");
  tree->Branch("min",min.data(), "min[nCh]/F");
  tree->Branch("t",t, "t[nCh]/F");
  tree->Branch("tSiPM", tSiPM, "tSiPM[nCh]/F");
  tree->Branch("BL", BL, "BL[nCh]/F");
  tree->Branch("BL_RMS", BL_RMS, "BL_RMS[nCh]/F");
  tree->Branch("Integral_0_300", Integral_0_300, "Integral_0_300[nCh]/F");
  tree->Branch("Integral", Integral, "Integral[nCh]/F");
  tree->Branch("Integral_Correction",&Integral_Correction, "Integral_Correction/F");
  tree->Branch("Integral_mVns", Integral_mVns, "Integral_mVns[nCh]/F");
  tree->Branch("EventIDsamIndex",EventIDsamIndex, "EventIDsamIndex[nCh]/I");
  tree->Branch("FirstCellToPlotsamIndex",FirstCellToPlotsamIndex, "FirstCellToPlotsamIndex[nCh]/I");

// // //   tree->Branch("MeasuredBaseline_usbwc", MeasuredBaseline_usbwc, baseline_ss.Data());
// // //   tree->Branch("AmplitudeValue_usbwc", AmplitudeValue_usbwc, amplitude_ss.Data());
// // //   tree->Branch("ComputedCharge_usbwc", ComputedCharge_usbwc, charge_ss.Data());
// // //   tree->Branch("RiseTimeInstant_usbwc", RiseTimeInstant_usbwc, leadingEdgeTime_ss.Data());
// // //   tree->Branch("FallTimeInstant_usbwc", FallTimeInstant_usbwc, trailingEdgeTime_ss.Data());
// // //   tree->Branch("RawTriggerRate_usbwc", RawTriggerRate_usbwc, rateCounter_ss.Data());
  //tree->Branch("amplValues", amplValues, "amplValues[nCh][1024]/S");
 // tree->Branch("hCh","TH1F",&hCh,128000,1);
  ///////////////////////////////////////////////////////

    int nitem = 1;
    ifstream inList;
    TString fileName;
    inList.open(_inFileList);
    assert(inList.is_open());

    int wavePrintStatus=-1;
    int ch0PrintStatus=-1;
    int trigPrintStatus=-1;
    int signalPrintStatus=-1;
    while(inList >> fileName){
      fileName = _inDataFolder + fileName;
      cout << endl;
      cout << fileName << endl;
      FILE* pFILE = fopen(fileName.Data(),"rb");
      if (pFILE==NULL) {fputs ("File error",stderr); assert(0);}
      //cout<<" ---> File to convert : " << fileName << endl;
      fseek (pFILE , 0 , SEEK_END);
      int totFileSizeByte = ftell (pFILE);
      rewind (pFILE);
      cout<<"totFileSizeByte = "<<totFileSizeByte<<endl;
      int size_of_header;
      if (WCVersion == WCHU){
        size_of_header = 328;
      }
      else if (WCVersion == WCAlexander){
        size_of_header = 327;
      }
      char header[size_of_header];
        nitem=fread(header,1,size_of_header,pFILE);
      // }
      // else if (WCVersion == WCAlexander){
      //   char header[327];
      //   nitem = fread(header,1,327,pFILE);
      // }
      cout << "Header:\n" << header << endl;

      char* word;
      word = strtok(header," \n");
      while(word != NULL){
	  if(strcmp("ACQUIRED:",word) == 0){
	    word = strtok(NULL, " \n");
	    nActiveCh = atoi(word);
	    break;
	  }
	 //printf ("%s\n",word);
	 word = strtok(NULL, " \n");
      }

      if(nActiveCh>9){
	cout << endl;
	char dummy;
	nitem=fread(&dummy,1,1,pFILE);
      }

      int whileCounter = 0;
      while(nitem>0){ //event loop
      std::vector<TObject*> eventTrash;
      
      whileCounter++;
	nitem = fread (&EventNumber	,sizeof(int), 1,pFILE);
	nitem = fread (&EpochTime	,sizeof(double)      , 1,pFILE);
	nitem = fread (&Year		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Month		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Day		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Hour		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Minute		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Second		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Millisecond	,sizeof(unsigned int), 1,pFILE);
  if (WCVersion == WCHU){
    nitem = fread (&nCh ,sizeof(unsigned int),1,pFILE); // since V2.8.14 the number of stored channels is written for each event
  }
  else if (WCVersion == WCAlexander){
    nCh = 16;
  }
	


	if(EventNumber%100==0)printf("POS, ev, y-m-d-h-min-s-ms, nActive-nCh: %ld, %d, %d-%d-%d-%d-%d-%d-%d, %d-%d \n", ftell(pFILE), EventNumber,Year,Month,Day,Hour,Minute,Second,Millisecond,nActiveCh,nCh);

	float	MeasuredBaseline[16];
	float	AmplitudeValue[16];
	float	ComputedCharge[16];
	float	RiseTimeInstant[16];
	float	FallTimeInstant[16];
	float	RawTriggerRate[16];
	float floatR=-1;
        for(int i = 0;i<nCh;i++){
	  //printf("i, currentPositionByte %d %ld\n",i,ftell(pFILE));
	  nitem = fread (&ChannelNr[i]	       ,sizeof(int),1,pFILE);
          nitem = fread (&EventIDsamIndex[i]        ,sizeof(int),1,pFILE);
	  nitem = fread (&FirstCellToPlotsamIndex[i],sizeof(int),1,pFILE);
	  nitem = fread (&floatR,1,4,pFILE); MeasuredBaseline[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); AmplitudeValue[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); ComputedCharge[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); RiseTimeInstant[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); FallTimeInstant[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); RawTriggerRate[i] = floatR;
	  ChannelNr[i]=i;

	  TString title("");
	  title.Form("ch %d, ev %d",i,EventNumber);
	  hCh.Reset();
	  hCh.SetTitle(title);
	  
	  if(i>=14){
	    for(int j = 0;j<1024;j++){
	      nitem = fread (&amplValues[i][j],sizeof(short),1,pFILE);
	      hCh.SetBinContent(j+1,-(amplValues[i][j]*coef*1000));
	    }//for 1024
	  }
	  else{
	    for(int j = 0;j<1024;j++){
	      nitem = fread (&amplValues[i][j],sizeof(short),1,pFILE);
	      hCh.SetBinContent(j+1,(amplValues[i][j]*coef*1000));
	    }//for 1024
	  }
	  
	  //for(int t=0;t<nActiveCh-nCh;t++){
	  //  int dummy;
	  //  nitem = fread(&dummy,sizeof(int),1,pFILE);
	  //  printf("trigger channel number: %d\n",dummy);
	  //}

	  max[i]=hCh.GetMaximum();
	  min[i]=hCh.GetMinimum();
	  amp[i]=hCh.GetMaximum();
	  
	  getBL(&hCh, BL_output,0,30);
	  BL[i] = BL_output[0];
	  BL_RMS[i] = BL_output[1];

    amp[i]=hCh.GetMaximum() - BL[i];
	  
	  for(int j=1;j<=hCh.GetXaxis()->GetNbins();j++){
	    hCh.SetBinError(j,BL_RMS[i]);
	  }
	  hChtemp.at(i) = hCh;
	  t[i]=CDF(&hCh,fTrigFit,0.1);
    if (i>6){
      t[i]=CDF(&hCh,fTrigFit,0.5);
    }

    if(EventNumber%wavesPrintRate==0){
      cWaves.cd(1+4*(i%4)+(i)/4);
      hCh.DrawCopy();
      TLine* ln = new TLine(t[i],-2000,t[i],2000);
      ln->SetLineColor(2);
      ln->Draw("same");
    }

	  if(i<=5||i==15)Integral_0_300[i] = (hCh.Integral(1, 1024, "width")-BL[i]*1024*SP)/pe;//Calculating Integral of histogram from 0 to 300ns; starting from bin 1 (0 is the overflow bin) to bin corresponding to 300ns. Option "width" multiplies by bin-width such that the integral is independant of the binning
	  else Integral_0_300[i] = (hCh.Integral(1, 1024, "width")-BL[i]*1024*SP);
	  
          if(EventNumber%ch0PrintRate==0&&i==0){
	    //cCh0.cd(1);
	    cCh0.cd();
	    hCh.DrawCopy("hist");
	    TLine ln(t[0],-2000,t[0],2000);
	    ln.SetLineColor(2);
	    //hCh.GetYaxis()->SetRangeUser(-10,1200);
	    ln.Draw("same");
	    if(ch0PrintStatus<0){cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf("),"pdf");ch0PrintStatus=0;}
	    else cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf"),"pdf");
	  }

	 if(EventNumber%trigPrintRate==0&&(i<4)){
	    cTrig.cd(i+1);
	    hCh.DrawCopy();
	    TLine* ln = new TLine(t[i],-2000,t[i],2000);
	    ln->SetLineColor(2);
	    ln->Draw("same");
	    //fTrigFit->DrawCopy("same");
	    eventTrash.push_back(ln);
	  }

	 if(EventNumber%signalPrintRate==0&&(i>=4&&i<=7)){
	    cSignal.cd(i+1-4);
	    hCh.DrawCopy();
	    TLine* ln = new TLine(t[i],-2000,t[i],2000); 
	    ln->SetLineColor(2);
	    ln->Draw("same");
	    eventTrash.push_back(ln);
	  }

       }//for nCh

      trigT = (t[0]+t[1]+t[2]+t[3])/4;
      trigTp = (t[0]+t[1]-t[2]-t[3])/4;
      t0t1 = (t[0]-t[1]);
      t2t3 = (t[2]-t[3]);
      tPMT1 = t[4]-trigT;
      tPMT2 = t[5]-trigT;
      for (int i=0; i<=7; i++){
        tSiPM[i+7] = t[i+7] - trigT;
        if (tSiPM[i+7] < -66){
          t[i+7] = CDFinvert(&hChtemp.at(i+7),0.33);
          tSiPM[i+7] = t[i+7] - trigT;
        }
      }
      if(tPMT2<-52){
        t[5]=CDFinvert(&hChtemp.at(5),0.1);
        tPMT2 = t[5]-trigT;
      }
      //tPMT2i = iCFD(&hChtemp.at(5),trigT-55,2,BL[5])-trigT;
      Integral[5] = integral(&hChtemp.at(5),t[5]-5,t[5]+65,BL[5])/pe;
      for (int i=0; i<=8; i++){
        Integral[i+6] = (integral(&hChtemp.at(i+6), trigT-65, trigT+5, BL[i+6]));
      }
      
      hChtemp.at(6).Add(&hChtemp.at(6),&hChtemp.at(7),1,1);
   //     if(EventNumber%wavesPrintRate==0){
	  //   cWaves.cd(1+4*(6%4)+(6)/4);
	  //   hChtemp.at(6).DrawCopy();
	  // }
      for (int i=0; i<=6; i++){
        Integral_mVns[i+8] = Integral[i+8] + SiPM_shift.at(i+1);
      }

        Integral_0_300[6] = Integral_0_300[6] + Integral_0_300[7];
      Integral_0_300[7] = Integral_0_300[6]-Integral_0_300[8]-Integral_0_300[9]-Integral_0_300[10]-Integral_0_300[11]-Integral_0_300[12]-Integral_0_300[13]-Integral_0_300[14];
      Integral[6] = Integral[6] + Integral[7];
      if (Integral[6]>14000){
      Integral_Correction = correction_function((Integral_mVns[8]+Integral_mVns[9]+Integral_mVns[10]+Integral_mVns[11]+Integral_mVns[12]+Integral_mVns[13]+Integral_mVns[14]));
    }
    else{
      Integral_Correction = 0;
    }  
      Integral[6] = Integral[6] + Integral_Correction;
      Integral[7] = Integral[6]-Integral[8]-Integral[9]-Integral[10]-Integral[11]-Integral[12]-Integral[13]-Integral[14];
      Integral_mVns[6] = Integral[6];
      //Conversion to Number of p.e.
      for (int i=0; i<=7; i++){
        Integral_0_300[i+7] = (Integral_0_300[i+7] + SiPM_shift.at(i)) / pe_SiPM.at(i);
        //Integral_mVns[i+7] = Integral[i+7] + SiPM_shift.at(i);
        Integral[i+7] = (Integral[i+7] + SiPM_shift.at(i)) / pe_SiPM.at(i);
      }
Integral_mVns[7] = Integral[7] + SiPM_shift.at(0); 

      hChtemp.at(7).Add(&hChtemp.at(6),&hChtemp.at(8),1,-1);
      hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(9),1,-1);
      hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(10),1,-1);
      hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(11),1,-1);
      hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(12),1,-1);
      hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(13),1,-1);
      hChtemp.at(7).Add(&hChtemp.at(7),&hChtemp.at(14),1,-1);
   //    if(EventNumber%wavesPrintRate==0){
	  //   cWaves.cd(1+4*(7%4)+(7)/4);
	  //   hChtemp.at(7).DrawCopy();
	  // }
      
      
      tSUMp = t[6] - trigT;
      tSUMm = t[7] - trigT;
      
      
      
      if(max[15]>5||min[15]<-5)isVeto=1;
      else isVeto = 0;
      trigGate = abs(*(std::max_element(t,t+4))-*(std::min_element(t,t+4)));  
      
      if(max[0]<1240&&max[1]<1240&&max[2]<1240&&max[3]<1240&&isVeto==0){
	isTrig=1;
	if(isTrig&&BL[0]<1.1&&BL[1]<1.1&&BL[2]<1.1&&BL[3]<1.1){
	  isTrig=1;
	  if(trigT<140&&trigT>90&&trigGate<10){
	    isTrig=1;
	  }
	  else isTrig=0;
	}
	else isTrig=0;
      }
      else isTrig=0;
      if(isTrig==1){
	int shift = (int)((140-trigT)/SP);
	for(int j=0;j<(int)hChtemp.size();j++){
	  hChSum.at(j)->Add(&hChtemp.at(j),1);
	  hChShift_temp.at(j).Reset();
	  for(int bin=1;bin<=hCh.GetXaxis()->GetNbins()-shift;bin++){
	   hChShift_temp.at(j).SetBinContent(shift+bin,hChtemp.at(j).GetBinContent(bin));
	  }
	  hChShift.at(j)->Add(&hChShift_temp.at(j),1);
	}
      }
      
      
      if(isTrig==1&&max[5]<1240)isGoodSignal_5=1;
      else isGoodSignal_5=0;
      
      if(EventNumber%wavesPrintRate==0){
	    //TString plotSaveName("");
	    //plotSaveName.Form("%s/wave-%d.png",plotSaveFolder.Data(),EventNumber);
	    if(wavePrintStatus<0){cWaves.Print((TString)(plotSaveFolder+"/waves.pdf("),"pdf");wavePrintStatus=0;}
	    else cWaves.Print((TString)(plotSaveFolder+"/waves.pdf"),"pdf");
      }
      if(EventNumber%trigPrintRate==0){
	    if(trigPrintStatus<0){cTrig.Print((TString)(plotSaveFolder+"/trig.pdf("),"pdf");trigPrintStatus=0;}
	    else cTrig.Print((TString)(plotSaveFolder+"/trig.pdf"),"pdf");
      }
      if(EventNumber%signalPrintRate==0){
	    if(signalPrintStatus<0){cSignal.Print((TString)(plotSaveFolder+"/signal.pdf("),"pdf");signalPrintStatus=0;}
	    else cSignal.Print((TString)(plotSaveFolder+"/signal.pdf"),"pdf");
      }

      tree->Fill();
      //cleanEventMemory(eventTrash);
      }//while events

    }//while
    inList.close();
    cWaves.Clear();
    cWaves.Print((TString)(plotSaveFolder+"/waves.pdf)"),"pdf");
    cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf)"),"pdf");
    cTrig.Print((TString)(plotSaveFolder+"/trig.pdf)"),"pdf");
    cSignal.Print((TString)(plotSaveFolder+"/signal.pdf)"),"pdf");

  rootFile = tree->GetCurrentFile();
  rootFile->Write();
  rootFile->Close();
}