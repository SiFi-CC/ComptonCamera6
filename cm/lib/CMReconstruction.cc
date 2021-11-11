#include "CLog.hh"

#include <TCanvas.h>
#include <TF2.h>
#include <TMath.h>
#include <TRandom.h>

#include "CmdLineConfig.hh"

#include "CMReconstruction.hh"
#include "CMSimulation.hh"
#include "DataStructConvert.hh"
#include "Smoothing.hh"
#include "Sources/PointSource.hh"

using SiFi::tools::convertHistogramToMatrix;
using SiFi::tools::convertMatrixToHistogram;
using SiFi::tools::unvectorizeMatrix;
using SiFi::tools::vectorizeMatrix;

CMReconstruction::~CMReconstruction() = default;

CMReconstruction::CMReconstruction(TString simulationFile) {
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

void CMReconstruction::FillHMatrix(int nGamma) {
  log->info("CMReconstruction::FillHMatrix");

  int nIterations = nGamma;
  for (int objBin = 0; objBin < fObjectCoords.NBins(); objBin++) {
    double done = (double)objBin/fObjectCoords.NBins();
    if(std::floor(1000*done) == 1000*done){
      log->info("{} %",done*100);
    }
    int objBinX, objBinY;
    std::tie(objBinX, objBinY) = fObjectCoords.BinXY(objBin);
    log->debug("Hmatrix ({}, {}) voxel", objBinX, objBinY);

    // TODO switch to angular bins
    Double_t x = fObject.GetXaxis()->GetBinCenter(objBinX + 1);
    Double_t y = fObject.GetYaxis()->GetBinCenter(objBinY + 1);
    // Double_t y = fObject.GetYaxis()->GetBinCenter(fObject.GetNbinsY() - objBinY);
    // log->info("Source coord ({}, {})", x, y);
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

  if (CmdLineOption::GetStringValue("Hmatrix")) {
    TString hfilename(CmdLineOption::GetStringValue("Hmatrix"));
    log->info("Hmatrix file: {}", hfilename);
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
    log->info("Mask1 bins: {}", fMaskCheckCoords.NBins());
    log->info("Mask2 bins: {}", fMaskCoords.NBins());

    Double_t ma1 = fMaskCheck.GetA(), mb1 = fMaskCheck.GetB(),
             mc1 = fMaskCheck.GetC();
    Double_t md1 = fMaskCheck.GetD(), mY1 = fMaskCheck.GetDimY(),
             mZ1 = fMaskCheck.GetDimZ();
    Double_t da1 = fDetPlaneCheck.GetA(), db1 = fDetPlaneCheck.GetB(),
             dc1 = fDetPlaneCheck.GetC();
    Double_t dd1 = fDetPlaneCheck.GetD(), dY1 = fDetPlaneCheck.GetDimY(),
             dZ1 = fDetPlaneCheck.GetDimZ();

    Double_t ma2 = fMask.GetA(), mb2 = fMask.GetB(), mc2 = fMask.GetC();
    Double_t md2 = fMask.GetD(), mY2 = fMask.GetDimY(), mZ2 = fMask.GetDimZ();
    Double_t da2 = fDetPlane.GetA(), db2 = fDetPlane.GetB(),
             dc2 = fDetPlane.GetC();
    Double_t dd2 = fDetPlane.GetD(), dY2 = fDetPlane.GetDimY(),
             dZ2 = fDetPlane.GetDimZ();

    if (fMaskCheckBins != fMaskBins) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Fmask HFile Nbins = {}", fMaskCheckBins);
      log->info("Fmask InputDataFile Nbins = {}", fMaskBins);
      exit(EXIT_FAILURE);
    } else if (ma1 != ma2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Fmask HFile A = {}", ma1);
      log->info("Fmask InputDataFile A = {}", ma2);
      exit(EXIT_FAILURE);
    } else if (mb1 != mb2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Fmask HFile B = {}", mb1);
      log->info("Fmask InputDataFile B = {}", mb2);
      exit(EXIT_FAILURE);
    } else if (mc1 != mc2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Fmask HFile C = {}", mc1);
      log->info("Fmask InputDataFile C = {}", mc2);
      exit(EXIT_FAILURE);
    } else if (md1 != md2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Fmask HFile D = {}", md1);
      log->info("Fmask InputDataFile D = {}", md2);
      exit(EXIT_FAILURE);
    } else if (mY1 != mY2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Fmask HFile DimY = {}", mY1);
      log->info("Fmask InputDataFile DimY = {}", mY2);
      exit(EXIT_FAILURE);
    } else if (mZ1 != mZ2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Fmask HFile DimZ = {}", mZ1);
      log->info("Fmask InputDataFile DimZ = {}", mZ2);
      exit(EXIT_FAILURE);
    }else if (da1 != da2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Detector HFile A = {}", da1);
      log->info("Detector InputDataFile A = {}", da2);
      exit(EXIT_FAILURE);
    } else if (db1 != db2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Detector HFile B = {}", db1);
      log->info("Detector InputDataFile B = {}", db2);
      exit(EXIT_FAILURE);
    } else if (dc1 != dc2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Detector HFile C = {}", dc1);
      log->info("Detector InputDataFile C = {}", dc2);
      exit(EXIT_FAILURE);
    } else if (dd1 != dd2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Detector HFile D = {}", dd1);
      log->info("Detector InputDataFile D = {}", dd2);
      exit(EXIT_FAILURE);
    } else if (dY1 != dY2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Detector HFile DimY = {}", dY1);
      log->info("Detector InputDataFile DimY = {}", dY2);
      exit(EXIT_FAILURE);
    } else if (dZ1 != dZ2) {
      log->error("Inconsistent parameters of H matrix and input data");
      log->info("Detector HFile DimZ = {}", dZ1);
      log->info("Detector InputDataFile DimZ = {}", dZ2);
      exit(EXIT_FAILURE);
    }

  } else {
    log->info("Hmatrix will be calculated with 100000 events");
    FillHMatrix(100000);
  }

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
  fMask.Write("mask");
  fDetPlane.Write("detector");
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

void CMReconstruction::HmatrixToFile(const TString& filename) {
  TFile file(filename, "RECREATE");
  file.cd();

  log->info("FILL Hmatrix");
  FillHMatrix(CmdLineOption::GetIntValue("Events"));

  log->info("WRITE Hmatrix");
  fMatrixH.Write("matrixH");
  fMask.Write("mask");
  fDetPlane.Write("detector");

  file.Close();
}
