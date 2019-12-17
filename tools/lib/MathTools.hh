#ifndef __MathTools_H
#define __MathTools_H 1

#include "TObject.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TF1.h"
#include "TH2.h"
#include "TString.h"
#include <iostream>

namespace MathTools{
  
    void ScaleGraphAndMove(TGraphErrors* g, Double_t factor=1, Double_t offset=0);
    void ScaleGraphAndMove(TGraph* g, Double_t factor=1, Double_t offset=0);
    TGraphErrors* SubtractBackground(TGraphErrors* gorig, TGraphErrors* bg);
    TGraphErrors* SubtractBackground(TGraphErrors* gorig, TF1* bg);
    TGraphErrors* SubtractBackground(TGraphErrors* gorig, TSpline3* bg);
    TH2 *SmoothGauss(TH2 *hin, double sigma);
    
};

#endif
