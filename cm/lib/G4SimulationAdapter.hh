#ifndef __G4SimulationAdapter_H_
#define __G4SimulationAdapter_H_ 1

#include "BinnedGeometry.hh"
#include "CLog.hh"
#include <TTree.h>
#include <vector>

struct CameraGeometry {
  BinnedGeometry source;
  BinnedGeometry mask;
  BinnedGeometry detector;
  std::vector<TFile*> recoData;

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

class G4SimulationAdapter : public TObject {
public:
  G4SimulationAdapter() = default;
  G4SimulationAdapter(TString filename);
  virtual ~G4SimulationAdapter();

  std::vector<TFile*> Filter(std::function<bool(TFile*)> filter);
  CameraGeometry GetFirstReconstructData();
  void VerifyForReconstruct(TFile* simulationFile);

private:
  void ReadMetadata();

  void ParseSelected(CameraGeometry* camera);
  void ParsePointSources(CameraGeometry* camera);
  bool IsSimulationGeometryEqual(TList* sim1, TList* sim2);

  std::vector<TFile*> fFiles = {};
  TFile* fSelected = nullptr;
  std::vector<TFile*> fSimulations;
  std::unordered_set<std::string> fIgnoredKeys = {"sourcePosX", "sourcePosY"};

  SiFi::logger log = SiFi::createLogger("G4SimulationAdapter");
  ClassDef(G4SimulationAdapter, 0)
};

#endif
