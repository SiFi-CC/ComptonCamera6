#include "CMReconstruction.hh"
#include "CLog.hh"
#include "CMSimulation.hh"
#include "DataStructConvert.hh"
#include "Sources/PointSource.hh"
#include <TCanvas.h>
#include <TF2.h>
#include <TMath.h>
#include <TRandom.h>
#include "CmdLineConfig.hh"

using SiFi::tools::convertHistogramToMatrix;
using SiFi::tools::convertMatrixToHistogram;
using SiFi::tools::unvectorizeMatrix;
using SiFi::tools::vectorizeMatrix;

CMReconstruction::~CMReconstruction() = default;

CMReconstruction::CMReconstruction(TString simulationFile) {
  // log->info("Load simulation results from file");
  log->info("Load simulation results from file {}" , simulationFile);
  TFile file(simulationFile, "READ");
  file.Print();

  log->debug("Load from \"mask\" object of type {}",
             file.Get("mask")->ClassName());
  fMask = *static_cast<Mask*>(file.Get("mask"));

  log->debug("Load from \"detector\" object of type {}",
             file.Get("detector")->ClassName());
  fDetPlane = *static_cast<DetPlane*>(file.Get("detector"));

  log->debug("Load from \"H2Detector\" object of type {}",
             file.Get("H2Detector")->ClassName());
  fImage = *static_cast<TH2F*>(file.Get("H2Detector"));
  fImage.SetDirectory(nullptr);

  log->debug("Load from \"H2Source\" object of type {}",
             file.Get("H2Source")->ClassName());
  fObject = *static_cast<TH2F*>(file.Get("H2Source"));
  fObject.SetDirectory(nullptr);
  log->debug("file closed");

  fObjectCoords = H2Coords(&fObject);
  fImageCoords = H2Coords(&fImage);
  log->debug("dimensions of reconstructed object [{}]",
             fObjectCoords.String().c_str());
  log->debug("dimensions of detector image [{}]",
             fImageCoords.String().c_str());

  fImageMat.ResizeTo(fImageCoords.NBins(), 1);
  fImageMat = vectorizeMatrix(convertHistogramToMatrix(&fImage));

  fMatrixH.ResizeTo(fImageCoords.NBins(), fObjectCoords.NBins());
  fMatrixHPrime.ResizeTo(fObjectCoords.NBins(), fImageCoords.NBins());
  fMatrixH = TMatrixT<Double_t>(fImageCoords.NBins(), fObjectCoords.NBins());
  fMatrixHPrime =
      TMatrixT<Double_t>(fObjectCoords.NBins(), fImageCoords.NBins());
}

void CMReconstruction::FillHMatrix() {
  log->info("CMReconstruction::FillHMatrix");

  int nIterations = 10000;
  for (int objBin = 0; objBin < fObjectCoords.NBins(); objBin++) {
    int objBinX, objBinY;
    std::tie(objBinX, objBinY) = fObjectCoords.BinXY(objBin);
    log->debug("Hmatrix ({}, {}) voxel", objBinX, objBinY);

    // TODO switch to angular bins
    Double_t x = fObject.GetXaxis()->GetBinCenter(objBinX + 1);
    Double_t y = fObject.GetYaxis()->GetBinCenter(objBinY + 1);

    {
      PointSource src(TVector3(0, y, x), 1);
      CMSimulation sim(&src, &fMask, &fDetPlane);
      sim.SetLogLevel(spdlog::level::warn);
      sim.RunSimulation(nIterations);
      auto imgVec = vectorizeMatrix(convertHistogramToMatrix(sim.GetImage()));
      for (int imgBin = 0; imgBin < fImageCoords.NBins(); imgBin++) {
        if (imgVec(imgBin, 0) == 0) {
          // fix for division by zero durring reconstruction
          fMatrixH(imgBin, objBin) = 1e-6 / nIterations;
        } else {
          fMatrixH(imgBin, objBin) = imgVec(imgBin, 0) / nIterations;
        }
      }
    }
  }
  fMatrixHPrime.Transpose(fMatrixH);
}

Bool_t CMReconstruction::CalculateS() {
  log->info("CMReconstruction::CalculateS");

  fNormS = std::vector<Double_t>(fObjectCoords.NBins());
  for (int i = 0; i < fObjectCoords.NBins(); i++) {
    fNormS[i] = 0;
    for (int j = 0; j < fImageCoords.NBins(); j++) {
      fNormS[i] += fMatrixH(j, i);
    }
  }

  return true;
}

void CMReconstruction::SingleIteration() {
  log->debug("CMReconstruction::SingleIteration()  iter={}",
             fRecoObject.size());

  // H * f_k
  // vector of correlations
  // i-th element of this vector contains information whether image resulting
  // from source postioned in i-th bin correlates to current iteration.
  auto hfProduct = fMatrixH * fRecoObject.back();
  log->debug("SingleIteration  H * f_k ({}, {})", hfProduct.GetNrows(),
             hfProduct.GetNcols());

  // Image / (H * f_k)
  // Image vector with inverse of correlations used as weights, the largest
  // weight will be atributed for image element that has the worst corelation
  // with current iteration.
  TMatrixT<Double_t> weightedImage(fImageCoords.NBins(), 1);
  for (int i = 0; i < fImageCoords.NBins(); i++) {
    weightedImage(i, 0) = fImageMat(i, 0) / hfProduct(i, 0);
  }
  log->debug("SingleIteration Image / (H * f_k) ({}, {})",
             weightedImage.GetNrows(), weightedImage.GetNcols());

  TMatrixT<Double_t> nextIteration = fMatrixHPrime * weightedImage;
  for (int i = 0; i < fObjectCoords.NBins(); i++) {
    nextIteration(i, 0) =
        nextIteration(i, 0) * fRecoObject.back()(i, 0) / fNormS[i];
  }

  fRecoObject.push_back(nextIteration);

  log->debug("end CMReconstruction::SingleIteration()  iter={}",
             fRecoObject.size() - 1);
}

void CMReconstruction::RunReconstruction(Int_t nIterations) {
  log->info("CMReconstruction::RunReconstruction({})", nIterations);
  if (nIterations > 100) {
    log->error("Too many iterations requested. Currently <100 feasible. \nFor "
               "more please adjust the code.");
    throw "too many iterations";
  }


  if(CmdLineOption::GetStringValue("Hmatrix")){ 
    TString hfilename(CmdLineOption::GetStringValue("Hmatrix"));
    log->info("Hmatrix file: {}",hfilename);
    TFile hfile(hfilename);
    hfile.cd();

    fMatrixH.Read("matrixH");
    fMatrixHPrime.Transpose(fMatrixH);
    Mask fMaskCheck = *static_cast<Mask*>(hfile.Get("mask"));
    DetPlane fDetPlaneCheck = *static_cast<DetPlane*>(hfile.Get("detector"));

    H2Coords fMaskCheckCoords(fMaskCheck.GetPattern());
    H2Coords fMaskCoords(fMask.GetPattern());
    Int_t fMaskCheckBins = fMaskCheckCoords.NBins();
    Int_t fMaskBins = fMaskCoords.NBins();
    log->info("Mask1 bins: {}",fMaskCheckCoords.NBins());
    log->info("Mask2 bins: {}",fMaskCoords.NBins());

    Double_t ma1 = fMaskCheck.GetA(), mb1 = fMaskCheck.GetB(), mc1 = fMaskCheck.GetC();
    Double_t md1 = fMaskCheck.GetD(), mY1 = fMaskCheck.GetDimY(), mZ1 = fMaskCheck.GetDimZ();
    Double_t da1 = fDetPlaneCheck.GetA(), db1 = fDetPlaneCheck.GetB(), dc1 = fDetPlaneCheck.GetC();
    Double_t dd1 = fDetPlaneCheck.GetD(), dY1 = fDetPlaneCheck.GetDimY(), dZ1 = fDetPlaneCheck.GetDimZ();

    Double_t ma2 = fMask.GetA(), mb2 = fMask.GetB(), mc2 = fMask.GetC();
    Double_t md2 = fMask.GetD(), mY2 = fMask.GetDimY(), mZ2 = fMask.GetDimZ();
    Double_t da2 = fDetPlane.GetA(), db2 = fDetPlane.GetB(), dc2 = fDetPlane.GetC();
    Double_t dd2 = fDetPlane.GetD(), dY2 = fDetPlane.GetDimY(), dZ2 = fDetPlane.GetDimZ();

    //TODO: separate conditions in order to recognize reason of error
    if (fMaskCheckBins != fMaskBins || ma1 != ma2 || mb1 != mb2 || mc1 != mc2 || md1 != md2 || mY1 != mY2 || mZ1 != mZ2){
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Fmask HFile: \n"
                " Nbins = {}, A = {}, B = {}, C = {}, D = {}, DimY = {}, DimZ = {} \n\n", 
                  fMaskCheckBins, ma1, mb1, mc1, md1, mY1, mZ1);
      log->info("Fdet HFile: \n"
                "A = {}, B = {}, C = {}, D = {}, DimY = {}, DimZ = {} \n\n", 
                  da1, db1, dc1, dd1, dY1, dZ1);
      log->info("Fmask InputDataFile: \n" 
                "NBins = {},  A = {}, B = {}, C = {}, D = {}, DimY = {}, DimZ = {} \n\n", 
                  fMaskBins, ma2, mb2, mc2, md2, mY2, mZ2);
      log->info("Fdet InputDataFile: \n"
                " A = {}, B = {}, C = {}, D = {}, DimY = {}, DimZ = {} \n\n",
                  da2, db2, dc2, dd2, dY2, dZ2);                        
      exit(EXIT_FAILURE);//exit(number = 'error code')
    }

  } else {
    log->info("Hmatrix will be calculated");
    FillHMatrix();
  }
  // return;

  CalculateS();

  fRecoObject.push_back(TMatrixT<Double_t>(fObjectCoords.NBins(), 1));
  fRecoObject[0] = 1;

  for (int i = 0; i < nIterations; i++) {
    SingleIteration();
  }
}

void CMReconstruction::Write(TString filename) const {
  log->info("CMReconstruction::Write({})", filename.Data());
  TFile file(filename, "RECREATE");
  file.cd();

  fObject.Write();
  fImage.Write();
  fMatrixH.Write("matrixH");
  SiFi::tools::convertMatrixToHistogram(
      "histH", "histogram of matrix H(probability matrix)", fMatrixH)
      .Write();

  int nIterations = fRecoObject.size() - 1;
  std::vector<TH2F*> histReco(nIterations);
  std::vector<TH1D*> histProjX(nIterations);
  std::vector<TH1D*> histProjY(nIterations);

  TCanvas can("MLEM2D", "MLEM2D", 1000, 1000);
  TCanvas canz("MLEM1DZ", "MLEM1DZ", 1000, 1000);
  TCanvas cany("MLEM1DY", "MLEM1DY", 1000, 1000);

  can.Divide(static_cast<int>(sqrt(nIterations)) + 1,
             static_cast<int>(sqrt(nIterations)) + 1);
  canz.Divide(static_cast<int>(sqrt(nIterations)) + 1,
              static_cast<int>(sqrt(nIterations) + 1));
  cany.Divide(static_cast<int>(sqrt(nIterations)) + 1,
              static_cast<int>(sqrt(nIterations)) + 1);

  log->debug("Save {} iterations", nIterations);
  for (int i = 0; i < nIterations; i++) {
    log->debug("saving iteration {}", i);

    auto recoIteration = convertMatrixToHistogram(
        "reco", TString::Format("iteration %d", i).Data(),
        unvectorizeMatrix(fRecoObject[i], fObjectCoords.NBinsX(),
                          fObjectCoords.NBinsY()));
    recoIteration.Write();

    // Histograms can be released after writing canvas to file so it's necessary
    // to store all references to release memory later.
    histReco[i] = static_cast<TH2F*>(recoIteration.Clone());
    histReco[i]->SetDirectory(nullptr);
    histProjX[i] = recoIteration.ProjectionX(
        TString::Format("projection Z, iteration %d", i).Data());
    histProjY[i] = recoIteration.ProjectionY(
        TString::Format("projection Y, iteration %d", i).Data());

    can.cd(i + 1);
    gPad->SetLogz(1);
    histReco[i]->Draw("colz");

    canz.cd(i + 1);
    histProjX[i]->Draw();

    cany.cd(i + 1);
    histProjY[i]->Draw();
  }

  log->debug("Write() canvas recnstruction iterations");
  can.Write();
  log->debug("Write() canvas Z projection");
  canz.Write();
  log->debug("Write() canvas Y projection");
  cany.Write();

  for (auto& object : histReco) {
    delete object;
  }
  for (auto& object : histProjX) {
    delete object;
  }
  for (auto& object : histProjY) {
    delete object;
  }

  file.Close();
  log->debug("end CMReconstruction::Write({})", filename.Data());
}


void CMReconstruction::HmatrixToFile(TString filename) {
    TFile file(filename, "RECREATE");
    file.cd();

    log->info("FILL Hmatrix");
    FillHMatrix();

    log->info("WRITE Hmatrix");
    fMatrixH.Write("matrixH");
    fMask.Write("mask");
    fDetPlane.Write("detector");

    file.Close();
}
