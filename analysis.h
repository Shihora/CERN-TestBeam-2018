#ifndef ANALYSIS
#define ANALYSIS

#include <TString.h>
#include <TPolyMarker.h>

using namespace std;

string checkFilename(TString filename);
float CDF(TH1F* hWave,float thr);
float CDFinvert(TH1F* hWave,float thr);
float Integrate_50ns(TH1F* hWave, float BL);
float integral(TH1F* hWave,float t1,float t2,float BL);
float* getBL(TH1F* hWave, float* BL, float t1, float t2);
float* getBL_fit(TH1F* hWave, float* BL_fit_output, float t1, float t2);
double correction_function(double x);
void peakfinder(TH1F *hWave, int nPeaks, int sigma, double thr, double *Xarray, double *Yarray, TPolyMarker *pfMarker, bool pfON);

#endif