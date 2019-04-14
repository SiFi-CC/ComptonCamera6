#ifndef __G4SimulationAdapter_H_
#define __G4SimulationAdapter_H_ 1

#include "CLog.hh"
#include <TTree.h>
#include <vector>

class G4SimulationAdapter : public TObject {
public:
  G4SimulationAdapter() = default;
  G4SimulationAdapter(TString filename);
  virtual ~G4SimulationAdapter();

  std::vector<TFile*> Filter(std::function<bool(TFile*)> filter);
  std::vector<TFile*> GetFirstReconstructData();

private:
  void ReadMetadata();
  TFile* fFile = nullptr;
  std::vector<TFile*> fSimulations;

  SiFi::logger log = SiFi::createLogger("G4SimulationAdapter");
  ClassDef(G4SimulationAdapter, 0)
};

#endif
