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
#include <TText.h>
#include <TSpectrum.h> // peakfinder
#include <TPolyMarker.h> // peakfinder
#include <TError.h> // root verbosity level
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

float SP = 0.3125; // ns per bin
float pe = 47.46;//mV*ns
vector<float> SiPM_shift = {2.679, 2.532, 3.594, 3.855, 3.354, 3.886, 3.865, 4.754};
vector<float> calib_amp = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; // dummy
vector<float> calib_int = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; // dummy
vector<float> const_BL = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // dummy
vector<float> const_BL_AB = {-2.05,-1.43,-1.39,-2.63,-2.42,-2.34,-1.36,-1.00,-3.10,-2.14,-1.76,-2.95,-0.73,-0.69,-0.99};
vector<float> const_BL_CD = {-0.08,-0.39,-1.47,-0.56,-0.59,-0.85,-1.44,-1.00,-3.45,-3.47,-0.99,-0.68,-0.35,-0.40,-0.55};
vector<float> calib_amp_AB = {6.748,6.16313,6.07082,6.68036,6.65783,6.37541,6.7711,6.85418,6.68469,6.58283,6.98329,6.97906,6.76493,6.75924,6.78279,1};
vector<float> calib_amp_CD = {4.738141,4.689474,4.553902,4.554155,4.545284,4.577300,4.746832,4.396243,4.217127, 4.344094,4.416440,4.678121,4.678319,4.633572,4.705655,1};
vector<float> calib_amp_AB_new = {6.225833, 5.681876, 5.520674, 5.982826, 6.179563, 6.041097, 6.416068, 6.072533, 5.697019, 5.452204, 5.798762, 6.023438, 5.794798, 5.796922, 5.869892,1}
vector<float> calib_amp_CD_new = {4.661835, 4.543407, 4.356386, 4.440221, 4.389425, 4.484233, 4.662783, 3.939226, 3.737358, 3.876855, 3.971836, 4.207425, 4.142826, 4.179445, 4.226949,1}
vector<float> calib_int_AB_new = {54.339372, 51.120311, 48.323768, 51.724367, 53.002368, 51.895161, 53.368556, 54.160940, 50.392792, 48.624219, 52.848405, 52.114772, 51.153844, 50.862783, 50.617176,1};
vector<float> calib_int_CD_new = {44.267965, 43.887981, 42.386506, 42.066467, 41.266592, 42.462270, 42.703238,37.792835, 36.457600, 37.816670, 37.611428, 39.822824, 39.078728, 39.895177, 39.592268,1};

int wavesPrintRate = 1000;
int sumWOMAPrintRate = 1000;
int sumWOMBPrintRate = 1000;
int ch0PrintRate = 1000000;
int trigPrintRate = 1000000;
int signalPrintRate = 100000;
double coef = 2.5 / (4096 * 10);
string WCHU ("AB"), WCAlexander ("CD");

//External Variables - mostly definded in main.C
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
  /*
  Function that is used by main(). read() does most of the analysis of the data. The raw .bin files 
  are read, converted into root-histograms and then the analysis calls are done. Most of the longer
  functions are defined in alanysis.C.
  The read() function the saves all the events (event by event) of that particular run to a root tree
  which is the saved in /runs/runName/out.root.
  */
  
  /*Create root-file and root-tree for data*/
  TFile *rootFile = new TFile( _outFile, "RECREATE");
  if (rootFile->IsZombie()) {
    cout << "PROBLEM with the initialization of the output ROOT ntuple "
    << _outFile << ": check that the path is correct!!!"
    << endl;
    exit(-1);
  }
  TTree *tree = new TTree("T", "USBWC Data Tree");
  TTree::SetBranchStyle(0);

  gStyle->SetLineScalePS(1); // export high resolution plots

  /*Declare & define the variables that are to be saved in the root-tree or that are used during the analysis.*/
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
  Float_t tsumWOMA_invCFD = -999;
  Float_t tsumWOMB_invCFD = -999;
  Float_t tsumWOMA_invCFD_wrtTrig = -999;
  Float_t tsumWOMB_invCFD_wrtTrig = -999;
  Int_t isLastEvt = -999;
  Int_t isGoodSignal_5 = -999;
  Float_t trigGate = -999;
  Int_t nCh = -1;
  int nActiveCh = -1;
  Int_t ChannelNr[16];
  Int_t WOMID[16];  //1=A, 2=B, 3=C, 4=D

  float PE_WOM1, PE_WOM2; // calibrated, baseline-shifted sum signal
  float PE_WOM1_int, PE_WOM2_int; // calibrated, baseline-shifted sum signal
  float t_PE_WOM1, t_PE_WOM2; // point in time of sum signal
  float chPE[16]; // single channel amplitude at sum signal
  float chPE_int[16]; // single channel integral

  std::vector<float> amp(16,-999);
  std::vector<float> amp_inRange(16,-999);
  std::vector<float> max(16,-999);
  std::vector<float> min(16,-999);
  Float_t t[16];
  Float_t tSiPM[16];

  float Integral_0_300[16];//array used to store Integral of signal from 0 to 300ns
  float Integral_inRange[16]; // calculate integral in given range
  float Integral[16];
  float Integral_mVns[16];

  float BL_output[4];//array used for output getBL-function
  Float_t BL_lower[16];//store baseline for 16 channels for 0-75ns range
  Float_t BL_RMS_lower[16];//store rms of baseline for 16 channels for 0-75ns range
  Float_t BL_Chi2_lower[16];//store chi2/dof of baseline-fit for 16 channels for 0-75ns range
  Float_t BL_pValue_lower[16];
  Float_t BL_upper[16];//store baseline for 16 channels for 220-320ns range
  Float_t BL_RMS_upper[16];//store rms of baseline for 16 channels for 220-320ns range
  Float_t BL_Chi2_upper[16];//store chi2/dof of baseline-fit for 16 channels for 220-320ns range
  Float_t BL_pValue_upper[16];

  Float_t BL_used[16];
  Float_t BL_Chi2_used[16];
  Float_t BL_pValue_used[16];
  float noiseLevel[16];

  int nPeaks = 4; // maximum number of peaks to be stored by peakfinder; has to be set also when creating branch
  Double_t peakX[16][nPeaks];
  Double_t peakY[16][nPeaks];

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
  // uncommtent, if .root file name should equal raw data file
  TString plotSaveFolder  = _inDataFolder;
  plotSaveFolder.ReplaceAll("data","runs");
  // TString plotSaveFolder  = _outFile;
  // plotSaveFolder.ReplaceAll("out.root","");
  TCanvas cWaves("cWaves","cWaves",1000,700);
  cWaves.Divide(4,4);
  TCanvas csumWOMA("csumWOMA","csumWOMA",1000,700);
  csumWOMA.Divide(4,2);
  TCanvas csumWOMB("csumWOMB","csumWOMB",1000,700);
  csumWOMB.Divide(3,3);
  TCanvas cCh0("cCh0","cCh0",1500,900);
  cCh0.Divide(2,2);
  TCanvas cTrig("cTrig","cTrig",1500,900);
  cTrig.Divide(2,2);
  TCanvas cSignal("cSignal","cSignal",1500,900);
  cSignal.Divide(2,2);
  TCanvas cChSum("cChSum","cChSum",1500,900);
  cChSum.Divide(4,4);

  /*Create branches in the root-tree for the data.*/
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
 
  // CHANNEL INFO (but everything that is nCH-dependend below)
  tree->Branch("nCh",&nCh, "nCh/I");
  tree->Branch("WOMID",WOMID,"WOMID[nCh]/I");
  tree->Branch("ch",ChannelNr, "ch[nCh]/I");
  // AMPLITUDE
  tree->Branch("amp",amp.data(), "amp[nCh]/F"); // calibrated
  tree->Branch("amp_inRange",amp_inRange.data(), "amp_inRange[nCh]/F"); // calibrated
  tree->Branch("max",max.data(), "max[nCh]/F");
  tree->Branch("min",min.data(), "min[nCh]/F");
  // INTEGRAL
  tree->Branch("Integral_0_300", Integral_0_300, "Integral_0_300[nCh]/F");
  tree->Branch("Integral_inRange", Integral_inRange, "Integral_inRange[nCh]/F");
  tree->Branch("Integral", Integral, "Integral[nCh]/F"); // calibrated
  tree->Branch("Integral_mVns", Integral_mVns, "Integral_mVns[nCh]/F"); // calibrated
  // TIMING
  tree->Branch("t",t, "t[nCh]/F");
  tree->Branch("tSiPM", tSiPM, "tSiPM[nCh]/F");
  // BASELINE
  tree->Branch("BL_lower", BL_lower, "BL_lower[nCh]/F");
  tree->Branch("BL_RMS_lower", BL_RMS_lower, "BL_RMS_lower[nCh]/F");
  tree->Branch("BL_Chi2_lower", BL_Chi2_lower, "BL_Chi2_lower[nCh]/F");
  tree->Branch("BL_pValue_lower", BL_pValue_lower, "BL_pValue_lower[nCh]/F");
  tree->Branch("BL_upper", BL_upper, "BL_upper[nCh]/F");
  tree->Branch("BL_RMS_upper", BL_RMS_upper, "BL_RMS_upper[nCh]/F");
  tree->Branch("BL_Chi2_upper", BL_Chi2_upper, "BL_Chi2_upper[nCh]/F");
  tree->Branch("BL_pValue_upper", BL_pValue_upper, "BL_pValue_upper[nCh]/F");
  tree->Branch("BL_used", BL_used, "BL_used[nCh]/F");
  tree->Branch("BL_Chi2_used", BL_Chi2_used, "BL_Chi2_used[nCh]/F");
  tree->Branch("BL_pValue_used", BL_pValue_used, "BL_pValue_used[nCh]/F");
  // PEAKFINDER
  tree->Branch("noiseLevel", noiseLevel, "noiseLevel[nCh]/F");
  tree->Branch("peakX",peakX,"peakX[nCh][4]/D");
  tree->Branch("peakY",peakY,"peakY[nCh][4]/D");
  // CALIBRATED SUM
  tree->Branch("PE_WOM1",&PE_WOM1, "PE_WOM1/F");
  tree->Branch("PE_WOM1_int",&PE_WOM1_int, "PE_WOM1_int/F");
  tree->Branch("PE_WOM2",&PE_WOM2, "PE_WOM2/F");
  tree->Branch("PE_WOM2_int",&PE_WOM2_int, "PE_WOM2_int/F");
  tree->Branch("t_PE_WOM1",&t_PE_WOM1, "t_PE_WOM1/F");
  tree->Branch("t_PE_WOM2",&t_PE_WOM2, "t_PE_WOM2/F");
  tree->Branch("chPE",chPE, "chPE[nCh]/F");
  tree->Branch("chPE_int",chPE_int, "chPE_int[nCh]/F");
  tree->Branch("tsumWOMA_invCFD",&tsumWOMA_invCFD,"tsumWOMA_invCFD/F");
  tree->Branch("tsumWOMB_invCFD",&tsumWOMB_invCFD,"tsumWOMB_invCFD/F");
  tree->Branch("tsumWOMA_invCFD_wrtTrig",&tsumWOMA_invCFD_wrtTrig,"tsumWOMA_invCFD_wrtTrig/F");
  tree->Branch("tsumWOMB_invCFD_wrtTrig",&tsumWOMB_invCFD_wrtTrig,"tsumWOMB_invCFD_wrtTrig/F"); 

  tree->Branch("EventIDsamIndex",EventIDsamIndex, "EventIDsamIndex[nCh]/I");
  tree->Branch("FirstCellToPlotsamIndex",FirstCellToPlotsamIndex, "FirstCellToPlotsamIndex[nCh]/I");

  /*Start reading the raw data from .bin files.*/
  int nitem = 1;
  ifstream inList;
  TString fileName;
  inList.open(_inFileList);
  assert(inList.is_open());

  int wavePrintStatus=-1;
  int sumWOMAPrintStatus=-1;
  int sumWOMBPrintStatus=-1;
  int ch0PrintStatus=-1;
  int trigPrintStatus=-1;
  int signalPrintStatus=-1;
  while(inList >> fileName){
    fileName = _inDataFolder + fileName;
    cout << endl;
    cout << fileName << endl;
    FILE* pFILE = fopen(fileName.Data(),"rb");
    if (pFILE==NULL) {
      fputs ("File error",stderr); 
      assert(0);
    }
    fseek (pFILE , 0 , SEEK_END);
    int totFileSizeByte = ftell (pFILE);
    rewind (pFILE);
    cout<<"totFileSizeByte = "<<totFileSizeByte<<endl;
    int size_of_header;
    /*During 2018 testbeam measurements two WaveCatchers were used. One from the Berlin group
    and one from Alexander from Geneva. As these two Wavecatchers had two different versions
    there are two types of raw data files that have different header lengths.*/
    if (WCVersion == WCHU){
      size_of_header = 328;
      calib_amp = calib_amp_AB_new;
      calib_int = calib_int_AB_new;
      const_BL = const_BL_AB;
    }
    else if (WCVersion == WCAlexander){
      size_of_header = 327;
      calib_amp = calib_amp_CD_new;
      calib_int = calib_int_CD_new;
      const_BL = const_BL_CD;
    }
    char header[size_of_header];
    nitem=fread(header,1,size_of_header,pFILE);

    cout << "Header:\n" << header << endl;

    char* word;
    word = strtok(header," \n");
    while(word != NULL){
      if(strcmp("ACQUIRED:",word) == 0){
        word = strtok(NULL, " \n");
        nActiveCh = atoi(word);
        break;
      }
      word = strtok(NULL, " \n");
    }

    if(nActiveCh>9){
      cout << endl;
      char dummy;
      nitem=fread(&dummy,1,1,pFILE);
    }

    int whileCounter = 0;
    /*Loop over events. Events are processed and analysed one by one in order.*/
    while(nitem>0){ //event loop
      std::vector<TObject*> eventTrash;
      whileCounter++;
      nitem = fread (&EventNumber, sizeof(int), 1, pFILE);
      nitem = fread (&EpochTime, sizeof(double), 1, pFILE);
      nitem = fread (&Year, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Month, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Day, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Hour, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Minute, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Second, sizeof(unsigned int), 1, pFILE);
      nitem = fread (&Millisecond, sizeof(unsigned int), 1, pFILE);
      if (WCVersion == WCHU){
        nitem = fread (&nCh ,sizeof(unsigned int),1,pFILE); // since V2.8.14 the number of stored channels is written for each event
      }
      else if (WCVersion == WCAlexander){
        nCh = 16;
      }
	


      if(EventNumber%100==0){
        printf("POS, ev, y-m-d-h-min-s-ms, nActive-nCh: %ld, %d, %d-%d-%d-%d-%d-%d-%d, %d-%d \n", ftell(pFILE), EventNumber,Year,Month,Day,Hour,Minute,Second,Millisecond,nActiveCh,nCh);
      }

      float	MeasuredBaseline[16];
      float	AmplitudeValue[16];
      float	ComputedCharge[16];
      float	RiseTimeInstant[16];
      float	FallTimeInstant[16];
      float	RawTriggerRate[16];
      float floatR=-1;

      /*Loop over individual channels. For each event the data from every channel is 
      processed and analysed one by one in order*/
      for(int i = 0;i<nCh;i++){
        nitem = fread (&ChannelNr[i], sizeof(int), 1, pFILE);
        nitem = fread (&EventIDsamIndex[i], sizeof(int), 1, pFILE);
        nitem = fread (&FirstCellToPlotsamIndex[i],sizeof(int),1,pFILE);
        nitem = fread (&floatR,1,4,pFILE); MeasuredBaseline[i] = floatR;
        nitem = fread (&floatR,1,4,pFILE); AmplitudeValue[i] = floatR;
        nitem = fread (&floatR,1,4,pFILE); ComputedCharge[i] = floatR;
        nitem = fread (&floatR,1,4,pFILE); RiseTimeInstant[i] = floatR;
        nitem = fread (&floatR,1,4,pFILE); FallTimeInstant[i] = floatR;
        nitem = fread (&floatR,1,4,pFILE); RawTriggerRate[i] = floatR;
        ChannelNr[i]=i;

        /*
        __ Set WOMID _________________________________________________________
        The labeling of the WOMs in the box was done using the letters A,B,C,D. For convinience these letters are here replaced by the numbers 1-4 which is stored in the root-tree for every channel and every event.
        */
        if (WCVersion == WCAlexander){
          if (i <= 6){ WOMID[i] = 3; }
          else if (i >= 7 && i <= 14){ WOMID[i] = 4; }
        }
        else {
          if (i <= 6){ WOMID[i] = 1; }
          else if (i >= 7 && i <= 14){ WOMID[i] = 2; }
        }

        TString title("");
        title.Form("ch %d, ev %d",i,EventNumber);
        hCh.Reset();
        hCh.SetTitle(title);

        /*
        __ Waveform Histogram _______________________________________________
        Writing the signal amplitude values into the root-histogram hCh.
        */
        if (i == 15){
          for(int j = 0;j<1024;j++){
            nitem = fread (&amplValues[i][j],sizeof(short),1,pFILE);
            hCh.SetBinContent(j+1,-(amplValues[i][j]*coef*1000));
          }
        }
        else {
          for(int j = 0;j<1024;j++){
            nitem = fread (&amplValues[i][j],sizeof(short),1,pFILE);
            hCh.SetBinContent(j+1,(amplValues[i][j]*coef*1000));
          }
        }

        /*The error of each value in each bin is set to 0.5 mV.*/
        for(int j=1;j<=hCh.GetXaxis()->GetNbins();j++){
          hCh.SetBinError(j,3);
        }

        /*Analysis if the event/signal starts.*/
        max[i] = hCh.GetMaximum();
        min[i] = hCh.GetMinimum();


        /*Saving the histogram of that event into a temporary histogram hChtemp. These histograms are available outside of the channel-loop. If analysis using the signals/events of multiple channels needs to be done, this can be accomplished by using hChtemp after the channel-loop.*/
        hChtemp.at(i) = hCh;
	  
        /*
        __ Baseline Fit _______________________________________________________
        Calculate baseline values infront and after the triggered signal
        Triggered signal is expected in the range fromm 100 to 150 ns
        */
        // BL_fit(&hChtemp.at(i), BL_output, 0.0, 75.0);
        BL_fit(&hChtemp.at(i), BL_output, 0.0, 30.0);
        BL_lower[i] = BL_output[0];
        BL_RMS_lower[i] = BL_output[1];
        BL_Chi2_lower[i] = BL_output[2];
        BL_pValue_lower[i] = BL_output[3];
        // BL_fit(&hChtemp.at(i), BL_output, 220.0, 320.0);
        BL_fit(&hChtemp.at(i), BL_output, 290.0, 320.0);
        BL_upper[i] = BL_output[0];
        BL_RMS_upper[i] = BL_output[1];
        BL_Chi2_upper[i] = BL_output[2];
        BL_pValue_upper[i] = BL_output[3];

        // determine "best" baseline
         if (BL_Chi2_upper[i] <= BL_Chi2_lower[i]){
          BL_used[i] = BL_upper[i];
          BL_Chi2_used[i] = BL_Chi2_upper[i];
          BL_pValue_used[i] = BL_pValue_upper[i];
        }
        else{
          BL_used[i] = BL_lower[i];
          BL_Chi2_used[i] = BL_Chi2_lower[i];
          BL_pValue_used[i] = BL_pValue_lower[i];
        }

        /*
        __ Peakfinder _________________________________________________________
        Implemented to search double-muon-event candiates
        Set maximum number of peaks stored in beginning of script -> nPeaks
        peakX/Yarray[nCh][nPeaks] stores peak coordinates as branches in tree
        Switch on/off with pfON
        -> when off:  set peakX/Yarray[nCh][nPeaks] to zero
        */
        gErrorIgnoreLevel = kError; // suppress root terminal output 

        bool pfON = false;
        if (i<15) {pfON = false;} // switch on/off peakfinder 
        int sigma = 10; // sigma of searched peaks
        Double_t thrPF = 0.1; // peakfinder threshold
        TPolyMarker pm; // store polymarker showing peak position, print later
        peakfinder(&hCh,0,130, nPeaks, sigma, thrPF, peakX[i], peakY[i], &pm, pfON);

        gErrorIgnoreLevel = kUnset; // return to normal terminal output

        // baseline-correct Y-values and convert to units of p.e.
        if (pfON)
        {
          for (int j = 0; j < nPeaks; ++j)
          {
            peakY[i][j] = amp2pe(peakY[i][j], calib_amp[i], BL_used[i]);
          }
        }

        // printf("X: %d %f %f %f %f \n",i,peakX[i][0],peakX[i][1],peakX[i][2],peakX[i][3]);
        // printf("Y: %d %f %f %f %f \n",i,peakY[i][0],peakY[i][1],peakY[i][2],peakY[i][3]);
        
        /*
        __ CFD _____________________________________________________________
        Setting the signal time by using a constant fraction disriminator method.
        The SiPM and the trigger sinals are handled differently using different thresholds.
        */
        if (i == 15){ //trigger
          t[i] = CDF(&hCh,0.5);
        }
        else { //SiPMs
          t[i] = CFD2(&hCh,0.35);
          if (t[i] < 95){
            t[i] = CFDinvert2(&hCh,0.35);
          }
        }

        /*
        __Print Raw Data to .txt ______________________________________________
        Select channel. Prints histogram raw data in a two column text file
        */

        // if (i==4 && BL_chi2[4]<1.7 && BL_chi2[4]>0.7)
        // if (i==4)
        // {
        //   TString histDataName;
        //   histDataName.Form("Ch%d_hist_data.txt",i);
        //   TString path2hist_data;
        //   path2hist_data.Form("%s/%s",(const char*)plotSaveFolder,(const char*)histDataName);
        //   FILE * histOut;
        //   histOut = fopen(path2hist_data,"a"); // produces overhead, maybe put this infront of loop
          
        //   Int_t nbins_x = hCh.GetNbinsX(); // bins at k==0 and k==nbins_x seem to have BinContent==0
        //   for (Int_t k=1; k<=nbins_x; k++)
        //   {
        //     fprintf(histOut,"%.4f %.8f\n",
        //     hCh.GetBinLowEdge(k)+hCh.GetBinWidth(k)/2,
        //     hCh.GetBinContent(k));
        //   }
        //   fclose(histOut);
        // }
        // clibrated sum
        if(EventNumber%sumWOMAPrintRate==0&&i<7){
            csumWOMA.cd(i+1);
            hCh.DrawCopy();
        }
        if(EventNumber%sumWOMBPrintRate==0&&i>6){
           csumWOMB.cd(i-6);
           hCh.DrawCopy();
        }

        /*
        __ Integral & Amplitude ________________________________________
        There are several definitions of the integral of a signal used here. Those are:
        - Integral_0_300: Integration over the entire time window (~320ns)
        - Integral: Integration over a smaller time window (~50ns) relative to the trigger
        Additionally the number of p.e. is now calculated using the amplitude
        and the calibration factors in the calib_amp-vactor. The function 'PE' calculates the amplitude of the signal, subtracts the better BL value and divides by the calibration factor.
        */

        Integral[i] = Integrate_50ns(&hCh, BL_used[i]) / calib_int.at(i); // difined 50 ns window
        Integral_inRange[i] = integral(&hCh, 100,125, BL_used[i]) / calib_int.at(i); // variable window

        // calibrated, BL-shifted amplitude at maximum in window
        amp[i] = PE(&hCh,calib_amp.at(i),BL_used[i], 100.0, 150.0);
        //maximum amplitude in range before expected signal (100-130 ns)
        amp_inRange[i] = PE(&hCh,calib_amp.at(i),BL_used[i], 0.0, 50.0);

        /*
        __ Printing Wafevorms ____________________________________________
        The signals for events can be printed to a .pdf file called waves.pdf. The rate at which the events are drawn to waves.pdf is set via the variable wavesPrintRate. Additional requirements can be set in the if-statement to look at specific events only.
        The entire if-statement so far also plots lines at the found signal maximum, the corresponding integration limit, as well as the BL values to each of the histograms.
        */
        if(EventNumber%wavesPrintRate==0){
          cWaves.cd(1+4*(i%4)+(i)/4);
          hCh.DrawCopy();
          hCh.GetXaxis()->SetRange((t[i]-20)/SP,(t[i]+30)/SP);
          int max_bin = hCh.GetMaximumBin();
          int lower_bin = max_bin - 20.0/SP;
          int upper_bin = max_bin + 30.0/SP;
          // double x = h->GetXaxis()->GetBinCenter(binmax);
          float max_time = hCh.GetXaxis()->GetBinCenter(max_bin);
          float lower_time = hCh.GetXaxis()->GetBinCenter(lower_bin);
          float upper_time = hCh.GetXaxis()->GetBinCenter(upper_bin);
          hCh.GetXaxis()->SetRange(0,1024);
          TLine* ln4 = new TLine(0,BL_lower[i],75,BL_lower[i]);
          TLine* ln5 = new TLine(220,BL_upper[i],320,BL_upper[i]);
          TText *text = new TText(.5,.5,Form("%f %f",BL_lower[i],BL_upper[i]));
          ln4->SetLineColor(2);
          ln5->SetLineColor(2);
          ln4->Draw("same");
          ln5->Draw("same");
          text->Draw("same");
          if (pfON){pm.Draw();} // print peakfinders polymarker
        }
      hCh.GetXaxis()->SetRange(1,30/SP);
      noiseLevel[i] = hCh.GetMaximum()-hCh.GetMinimum();
      hCh.GetXaxis()->SetRange(1,1024);
      // End of loop over inividual channels
      }

      /*
      __ Number of P.E. _____________________________________________________
      Calculate & save the number of p.e. for an entire WOM.
      Note: for the 1st WOM of each WC one channel was not recorded.
      Thus, there are only 7 values from 8. The result for the WOM is therefore
      scaled up by 8/7 to make the numbers comparable.
      */
      PE_WOM1 = 8/7*(amp[0]+amp[1]+amp[2]+amp[3]+amp[4]+amp[5]+amp[6]);
      PE_WOM2 = (amp[7]+amp[8]+amp[9]+amp[10]+amp[11]+amp[12]+amp[13]+amp[14]);

      /*
      __ TIMING _____
      */
      trigT = t[15];
      for (int i=0; i<=14; i++){
        tSiPM[i] = t[i] - trigT;
        /*
        if (tSiPM[i+7] < -66){
          t[i+7] = CDFinvert(&hChtemp.at(i+7),0.33);
          tSiPM[i+7] = t[i+7] - trigT;
        }
        */
      }

      /*
      __ FILLING SUM HISTOGRAMS ________________________
      For WOM 1 and 2 and determine time Resolution
      */

      // SUM WOM 1
      TH1F hSumA("hSumA","Sum A;ns;Amplitude, N_pe",1024,-0.5*SP,1023.5*SP);
      for(int hSumIndexA=0;hSumIndexA<7;hSumIndexA++)
      {
        TF1* f_const = new TF1("f_const","pol0",0,320);
        f_const->SetParameter(0,const_BL[hSumIndexA]);

        hChtemp.at(hSumIndexA).Add(f_const, -1);
        hChtemp.at(hSumIndexA).Scale(1.0/calib_amp.at(hSumIndexA));

        hSumA.Add(&hChtemp.at(hSumIndexA),1);
      }

      // get point of amplitude maximum
      PE_WOM1 = max_inRange(&hSumA, 100.0, 150.0)*8/7;
      t_PE_WOM1 = t_max_inRange(&hSumA, 100.0, 150.0);

      csumWOMA.cd(8);
      hSumA.DrawCopy();

      // get single channel amplitude and integral at time of sum maximum
      for (int i=0;i<7;i++)
      {
        chPE[i] = amp_atTime(&hChtemp.at(i), t_PE_WOM1);
        // reverse amplitude calibration before integration
        hChtemp.at(i).Scale(calib_amp.at(i));
        chPE_int[i] = integral(&hChtemp.at(i), t_PE_WOM1-10, t_PE_WOM1+15, 0)/calib_int.at(i);
      }

      PE_WOM1_int = chPE_int[0]+chPE_int[1]+chPE_int[2]+chPE_int[3]+chPE_int[4]+chPE_int[5]+chPE_int[6];

      // timing WOM 1
      tsumWOMA_invCFD = CFDinvert2(&hSumA,0.4);
      tsumWOMA_invCFD_wrtTrig = trigT-tsumWOMA_invCFD;

      // SUM WOM 2
      TH1F hSumB("hSumB","Sum B;ns;Amplitude, N_pe",1024,-0.5*SP,1023.5*SP);
      for(int hSumIndexB=7;hSumIndexB<15;hSumIndexB++)
      {
        TF1* f_const = new TF1("f_const","pol0",0,320);
        f_const->SetParameter(0,const_BL[hSumIndexB]);

        hChtemp.at(hSumIndexB).Add(f_const, -1);
        hChtemp.at(hSumIndexB).Scale(1.0/calib_amp.at(hSumIndexB));

        hSumB.Add(&hChtemp.at(hSumIndexB),1);
      }

      // get point of amplitude maximum
      PE_WOM2 = max_inRange(&hSumB, 100.0, 150.0);
      t_PE_WOM2 = t_max_inRange(&hSumB, 100.0, 150.0);

      csumWOMB.cd(9);
      hSumB.DrawCopy();

      // get single channel amplitude and integral at time of sum maximum
      for (int i=7;i<15;i++)
      {
        chPE[i] = amp_atTime(&hChtemp.at(i), t_PE_WOM2);
        // reverse amplitude calibration before integration
        hChtemp.at(i).Scale(calib_amp.at(i));
        chPE_int[i] = integral(&hChtemp.at(i), t_PE_WOM2-10, t_PE_WOM2+15, 0)/calib_int.at(i);
      }
      PE_WOM2_int = chPE_int[7]+chPE_int[8]+chPE_int[9]+chPE_int[10]+chPE_int[11]+chPE_int[12]+chPE_int[13]+chPE_int[14];

      // timing WOM 2
      tsumWOMB_invCFD = CFDinvert2(&hSumB,0.4);
      tsumWOMB_invCFD_wrtTrig = trigT-tsumWOMB_invCFD;
      
      // add-up all events channel-wise, not calibrated
      for (int i=0;i<=15;i++){
      	hChSum.at(i)->Add(&hChtemp.at(i),1);
      }

      /*
      */
      /* end */

      /*temporary for cfdScan_
      if (PE_WOM2 >= 5){    //Switch PE_WOM of you want WOM 1 or 2
        for(int i=0;i<N_CFD_points;i++){
          v_hTimeRes.at(i)->Fill(trigT-CFD2(&hSumB,0.01*i)); // Switch hSumA or B or CFD2 or CFDinvert2
        }
      }
      *///temporary end
      
      /*
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
      */

      /*Saving the plotted signals/events to a new page in the .pdf file.*/
      if(EventNumber%wavesPrintRate==0) {
        if(wavePrintStatus<0){
          cWaves.Print((TString)(plotSaveFolder+"/waves.pdf("),"pdf");
          wavePrintStatus=0;
        }
        else cWaves.Print((TString)(plotSaveFolder+"/waves.pdf"),"pdf");
      }
      if(EventNumber%sumWOMAPrintRate==0){
        if(sumWOMAPrintStatus<0){
          csumWOMA.Print((TString)(plotSaveFolder+"/sumWOMA.pdf("),"pdf");
          sumWOMAPrintStatus=0;
        }
        else csumWOMA.Print((TString)(plotSaveFolder+"/sumWOMA.pdf"),"pdf");
      }
      if(EventNumber%sumWOMBPrintRate==0){
        if(sumWOMBPrintStatus<0){
           csumWOMB.Print((TString)(plotSaveFolder+"/sumWOMB.pdf("),"pdf");
           sumWOMBPrintStatus=0;
        }
        else csumWOMB.Print((TString)(plotSaveFolder+"/sumWOMB.pdf"),"pdf");
      }
      if(EventNumber%trigPrintRate==0){
        if(trigPrintStatus<0){
          cTrig.Print((TString)(plotSaveFolder+"/trig.pdf("),"pdf");
          trigPrintStatus=0;
        }
        else cTrig.Print((TString)(plotSaveFolder+"/trig.pdf"),"pdf");
      }
      if(EventNumber%signalPrintRate==0){
        if(signalPrintStatus<0){
          cSignal.Print((TString)(plotSaveFolder+"/signal.pdf("),"pdf");
          signalPrintStatus=0;
        }
        else cSignal.Print((TString)(plotSaveFolder+"/signal.pdf"),"pdf");
      }

      /*Writing the data for that event to the tree.*/
      tree->Fill();
    }
    fclose(pFILE);
  }

  /*Clearing objects and saving files.*/
  inList.close();
  cWaves.Clear();
  cWaves.Print((TString)(plotSaveFolder+"/waves.pdf)"),"pdf");
  csumWOMA.Clear();
  csumWOMA.Print((TString)(plotSaveFolder+"/sumWOMA.pdf)"),"pdf");
  csumWOMB.Clear();
  csumWOMB.Print((TString)(plotSaveFolder+"/sumWOMB.pdf)"),"pdf");
  cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf)"),"pdf");
  cTrig.Print((TString)(plotSaveFolder+"/trig.pdf)"),"pdf");
  cSignal.Print((TString)(plotSaveFolder+"/signal.pdf)"),"pdf");
  for (int i=0; i<=15; i++){
  	cChSum.cd(i+1);
  	hChSum.at(i)->Draw();
  }
  cChSum.Print((TString)(plotSaveFolder+"/ChSum.pdf"),"pdf");

  /* temporary for cfdScan
  TH1F* sigma_timeRes_Fit = new TH1F("sigma_timeRes_Fit",";CFD threshold;Timeresolution, ns",N_CFD_points,0,0.01*N_CFD_points);
  TF1* fGaus = new TF1("fGaus","gaus",45,60);
  fGaus->SetLineWidth(1);
  TCanvas c_TimeRes("c_TimeRes");
  for(int i=0;i<N_CFD_points;i++){
    c_TimeRes.cd();
    v_hTimeRes.at(i)->Draw();
    v_hTimeRes.at(i)->Fit("fGaus","R");
    sigma_timeRes_Fit->SetBinContent(i+1,fGaus->GetParameter(2));
    sigma_timeRes_Fit->SetBinError(i+1,fGaus->GetParError(2));
    TString name("");
    name.Form("hTimeRes%d.png",i);
    c_TimeRes.Print((TString)(plotSaveFolder+"/"+name));
  }
  c_TimeRes.cd();
  sigma_timeRes_Fit->SetMarkerSize(1);
  sigma_timeRes_Fit->SetMarkerStyle(20);
  sigma_timeRes_Fit->SetLineColor(kRed);
  //c_TimeRes.Print((TString)("/home/maximilian/Dokumente/Arbeit/CERNTestBeam2018/Auswertung/CERN-TestBeam-2018/cfd-scan/"+name));
  c_TimeRes.Write();
  */// temporary end

  rootFile = tree->GetCurrentFile();
  rootFile->Write();
  rootFile->Close();
}