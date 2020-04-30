#ifndef __G4SimulationAdapter_H_
#define __G4SimulationAdapter_H_ 1

#include "BinnedGeometry.hh"
#include "CLog.hh"
#include <TTree.h>
#include <vector>

// Represents information about geometry of entire setup extrated
// from set of simulations.
struct CameraGeometry {
  BinnedGeometry source;
  BinnedGeometry mask;
  BinnedGeometry detector;
  // list of data files containg results of simulation for single point source
  // each.
  std::vector<TFile*> recoData;

  TMatrixT<Double_t> fMatrixHCam;

  void Print() {
    spdlog::info("Reconstructing for geometry");
    spdlog::info("Source");
    source.Print();
    spdlog::info("Mask");
    mask.Print();
    spdlog::info("Detector");
    detector.Print();
  }
};

// G4SimulationAdapter abstraction layer between raw root files and
// reconstruction alogrithm.
class G4SimulationAdapter : public TObject {
public:
  G4SimulationAdapter() = default;
  // Open file passed as arguemnt and all named filename.1, filename.2 ... that
  // exists. All simulations data present in those files are merged into single
  // list
  G4SimulationAdapter(TString filename);
  virtual ~G4SimulationAdapter();

  // Get list of simulation filtered with some condition
  std::vector<TFile*> Filter(std::function<bool(TFile*)> filter);
  // Gather data required for reconstuction, assume that first file defines
  // geometry of the system.
  CameraGeometry GetFirstReconstructData();
  // Check if reconsted simulaion had the same geomtry as reconstruct data.
  void VerifyForReconstruct(TFile* simulationFile);

private:
  // void ReadMetadata();

  // extract gemetry information from single simulation data file (fSelected)
  // and write it to camera object
  void ParseSelected(CameraGeometry* camera);
  // extract source geomtry by checking positions of all point sources and
  // assuming grid layout of those points
  // void ParsePointSources(CameraGeometry* camera);
  // verifies wheteher two simulations have the same geometry
  bool IsSimulationGeometryEqual(TList* sim1, TList* sim2);

  // list of data files
  std::vector<TFile*> fFiles = {};
  TFile* fSelected = nullptr;
  // list of signle point source simulations
  std::vector<TFile*> fSimulations;
  // metadata keys ignored when checking if two simulations have the same setup
  std::unordered_set<std::string> fIgnoredKeys = {"sourcePosX", "sourcePosY"};

  SiFi::logger log = SiFi::createLogger("G4SimulationAdapter");
  ClassDef(G4SimulationAdapter, 0)
};

#endif
