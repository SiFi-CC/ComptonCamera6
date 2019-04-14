#include "G4SimulationAdapter.hh"
#include <TFile.h>
#include <TParameter.h>

G4SimulationAdapter::G4SimulationAdapter(TString filename) {
  fFile = new TFile(filename, "READ");
  ReadMetadata();
}

G4SimulationAdapter::~G4SimulationAdapter() { delete fFile; }

void G4SimulationAdapter::ReadMetadata() {
  if (fFile->IsFolder()) {
    log->info("Found {} data files", fFile->GetListOfKeys()->GetSize());
  }
  auto fileContent = fFile->GetListOfKeys();
  for (auto element : *fileContent) {
    auto el = static_cast<TKey*>(element);
    fSimulations.push_back(static_cast<TFile*>(fFile->Get(el->GetName())));
  }
}

std::vector<TFile*>
G4SimulationAdapter::Filter(std::function<bool(TFile*)> filter) {
  std::vector<TFile*> filtered;
  std::copy_if(fSimulations.begin(), fSimulations.end(),
               std::back_inserter(filtered), filter);
  return filtered;
}

/**
 *  Potentially file could contain data for multiple simulations, e.g. different
 * energies, distance to mask, distance to detector, material etc. This function
 * will pick whathever is first on the list and find all matching to it.
 */
std::vector<TFile*> G4SimulationAdapter::GetFirstReconstructData() {
  auto first = *fSimulations.begin();
  auto firstUserInfo =
      static_cast<TTree*>(first->Get("metadata"))->GetUserInfo();

  log->info("Picked reconstruction data for:");
  for (auto entryFromFirst : *firstUserInfo) {
    // TODO: need to support string values in metadata
    if (TString(entryFromFirst->ClassName()) != "TParameter<double>") {
      continue;
    }
    auto e = static_cast<TParameter<double>*>(entryFromFirst);
    if (TString(e->GetName()) == "sourcePosX" ||
        TString(e->GetName()) == "sourcePosY") {
      continue;
    }
    log->info("key={}, value={}", e->GetName(), e->GetVal());
  }
  return Filter([firstUserInfo, this](TFile* file) -> bool {
    /* very inefficient implementation */
    auto userInfo = static_cast<TTree*>(file->Get("metadata"))->GetUserInfo();

    for (auto entryFromFirst : *firstUserInfo) {
      if (TString(entryFromFirst->ClassName()) != "TParameter<double>") {
        continue;
      }
      auto e1 = static_cast<TParameter<double>*>(entryFromFirst);
      if (TString(e1->GetName()) == "sourcePosX" ||
          TString(e1->GetName()) == "sourcePosY") {
        continue;
      }

      for (auto entry : *userInfo) {
        if (TString(entry->ClassName()) != "TParameter<double>") { continue; }
        auto e2 = static_cast<TParameter<double>*>(entry);
        if (TString(e1->GetName()) == TString(e2->GetName()) &&
            e1->GetVal() != e2->GetVal()) {
          log->debug("skipping because key={} have values ({}, {})",
                     e1->GetName(), e1->GetVal(), e2->GetVal());
          return false;
        }
      }
    }
    return true;
  });
}

