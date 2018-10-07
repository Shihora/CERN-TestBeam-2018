//root
#include <TString.h>
#include <TH1F.h>

//C, C++
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <iostream>

//local
#include "geometry.h"
#include "analysis.h"
#include "read.h"

using namespace std;

/* Declarations of global/external variables which are used accross all 
other files (read.C, geometry.C, analysis.C) */
string WCVersion;
int runNr = -999;
int pdgID= -999;
int isSP = -999;
int mp = -999;
int safPMT2 = -999;//solid angle factor of pmt 2
int safPMT1 = -999;//solid angle factor of pmt 1
int safSiPM = -999;//solid angle factor of SiPM
int trackL = -999;//track length
float energy = -999; // [GeV]
float horizontal= -999;
float vertical= -999;
float angle = -999;
std::vector<float> pmt2Pos = {410,410};


int main(int argc, char *argv[]){
  TString inFileList;
  TString inDataFolder;
  TString outFile;


  if(argc == 5){

    /* Used i.e. for calibration measurements which do not include angluar or 
    positional information as in testbeam measurements. Only information 
    on the local directories is handet to main() */

    inFileList = argv[1];
    inDataFolder = argv[2];
    outFile = argv[3];
    runNr=atoi(argv[4]);
    WCVersion = checkFilename(outFile);
    cout << WCVersion << endl;
    // if (WCVersion == 1){
    //   cout << "AB" << endl;
    //   return 0;
    // }
    // else if (WCVersion == 2){
    //   cout << "CD" << endl;
    //   return 0;
    // }
    cout<<"In data file list : "<<inFileList<<endl
      <<"In data path      : "<<inDataFolder<<endl
      <<"Out root file     : "<<outFile<<endl;
    read(inFileList, inDataFolder, outFile);
  }
  else if(argc == 13){

    /* Used for CERNTestBeam2017 data. The data from the testbeam has more 
    parameters which are handed to main() in order to calculate geometries 
    and save measurement positions, angles ect. */
    
    inFileList = argv[1];
    inDataFolder = argv[2];
    outFile = argv[3];
    runNr=atoi(argv[4]);
    horizontal = atof(argv[5])/1000; //units: [cm]
    vertical = atof(argv[6])/1000; //units: [cm]
    angle = atof(argv[7]);
    //argv[8] - nevents from table
    pdgID=atoi(argv[9]);// - pdgID
    energy=atof(argv[10]);// - energy [GeV]
    isSP=atoi(argv[11]);// - isSteelPlate
    mp=atoi(argv[12]);//  measure-point


    cout<<"In data file list : "<<inFileList<<endl
      <<"In data path      : "<<inDataFolder<<endl
      <<"Out root file     : "<<outFile<<endl
      <<"Run number         : "<<runNr<<endl;
    printf("hor,ver,ang: %4.2f %4.2f %4.2f\n",horizontal,vertical,angle);
    std::vector<float> starPos = getStartPos(horizontal,vertical,angle);
    //printf("start pos: %4.2f %4.2f %4.2f %4.2f\n",starPos[0],starPos[1],starPos[2],starPos[3]);
    std::vector<float> saf = solidAngleFactor(starPos,pmt2Pos);
    //printf("length: %4.2f %4.2f \n", saf[1],10*calculateDistance(horizontal,angle*TMath::Pi()/180));
    printf("length: %4.2f \n", saf[1]);
    printf("saf: %4.2f \n", saf[0]);
    safPMT2 = saf[0];
    trackL = saf[1];
	
    read(inFileList, inDataFolder, outFile);

  }
  else{
    cout<<" ERROR --->  in input arguments "<<endl
      <<"        [1] - in data file list"<<endl
      <<"        [2] - in data path"<<endl
      <<"        [3] - out root file"<<endl;
  }
  return 0;
}