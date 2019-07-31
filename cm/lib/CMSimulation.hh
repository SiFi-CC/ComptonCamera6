#ifndef __CMSimulation_H_
#define __CMSimulation_H_ 1
#include "CLog.hh"
#include "Source.hh"
#include "Track.hh"
#include <Mask.hh>
#include <TFile.h>
#include <TGeoManager.h>
#include <TH2F.h>
#include <TMatrix.h>
#include <TTree.h>

class CMSimulation : public TObject {

public:
  CMSimulation() = default;
  CMSimulation(Source* source, Mask* mask, DetPlane* detector)
      : fSource(source), fMask(mask), fDetPlane(detector) {
    Init();
  }
  virtual ~CMSimulation();
  void RunSimulation(Int_t nEvents);
  void ResetSimulation();

  void Write(TString name) const;
  void Print() const;

  TH2F* GetObject() const { return fH2Source; };
  TH2F* GetImage() const { return fH2Detector; };

  void SetLogLevel(spdlog::level::level_enum level) { log->set_level(level); }

private:
  void Init();
  Bool_t ProcessEvent();
  void BuildTGeometry(TString name) const;

  /** Source of radiation used in simulation */
  Source* fSource = nullptr;
  /** Mask implementation */
  Mask* fMask = nullptr;
  /** Detector plane */
  DetPlane* fDetPlane = nullptr;

  /** Stores data of simulated particles. Particle is saved only if projected
   * path is crossing mask and detector, otherwise is not scored anywhere.
   */
  TTree* fTree;
  /** Histogram representing distribution of tracks in Source
   *
   * TODO: This data will be moved to Source in the future.
   */
  TH2F* fH2Source;
  /** Histogram representing distribution of tracks in detector. This image is
   * used in reconstruction
   *
   * TODO: This data will be moved to Detector in the future.
   */
  TH2F* fH2Detector;
  /** Distribution of theta angles of registered particles
   *
   * TODO: This data will be removed or moved to Source in the future.
   */
  TH1F* fH1Theta;
  /** Distribution of theta angles of all generated particles
   *
   * TODO: This data will be removed or moved to Source in the future.
   */
  TH1F* fH1ThetaAll;
  /** Distribution of phi angles of registered particles
   *
   * TODO: This data will be removed or moved to Source in the future.
   */
  TH1F* fH1Phi;
  /** Distribution of phi angles of all generated particles
   *
   * TODO: This data will be removed or moved to Source in the future.
   */
  TH1F* fH1PhiAll;

  /**
   * Information about currently processed particle, needs to be kept in
   * object because TTree api requires source of data to be defined on init.
   * fTree object is Filled usng fields of this struct.
   */
  struct {
    /** Tracks representing paricle state at source mask an detector plane */
    Track sourceTrack, maskTrack, detectorTrack;
    /** Flag that describes whether particle was absorbed by mask */
    Bool_t absorbed;
  } fPersist;

  SiFi::logger log = SiFi::createLogger("CMSimulation");

  ClassDef(CMSimulation, 0)
};

#endif
