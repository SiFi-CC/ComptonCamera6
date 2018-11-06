#ifndef G4INPUT_H
#define G4INPUT_H

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#include <utility>
#include <iostream>
#include <sstream>
#include <stdlib.h> 
#include <vector>

#include "DR_GenerallStructs.hh"

/*! G4Input class. The class loads the output file of the G4Simulation and provides the different information:
 * - The information about the simulated CC dimensions
 * - The information about the real events
 * - The information about the reconstructed events.
 *One instance is build per run, the events can be accessed by LoadEvent() or by LoadEvent(Int_t no). */

class G4Input{
  public:
	/*! Constructer  
 * 	@param  InputFile - name of the file that shall be loaded.
 * 	@param  verb - set to true to get cout stream information.*/
	G4Input(TString InputFile,Bool_t verb);
	/*! Default Deconstructor*/
	~G4Input();
	/**@return fScatterer_x - The complete length of the scatterer in x dimension*/
	Double_t GetScattererXLength();
	/**@return fScatterer_y - The complete length of the scatterer in y dimension*/
	Double_t GetScattererYLength();
	/**@return fScatterer_z - The complete length of the scatterer in z dimension*/
	Double_t GetScattererZLength();
	/**@return fAbsorber_x - The complete length of the absorber in x dimension*/
	Double_t GetAbsorberXLength();
	/**@return fAbsorber_x - The complete length of the absorber in x dimension*/
	Double_t GetAbsorberYLength();
	/**@return fAbsorber_x - The complete length of the absorber in x dimension*/
	Double_t GetAbsorberZLength();
	/**@return fPositions_scatterer - The position of the scatterer*/
	TVector3 GetScattererPosition();
	/**@return fPositions_absorber - The position of the absorber*/
	TVector3 GetAbsorberPosition();
	/**@return the amount of Events that actually contain a compton event in the detector*/
	Int_t GetNumberOfRealEvents();
	/**@return the amount of Events that are stored in the input file*/
	Int_t GetNumberOfEvents();
	/**@return the amount of Events that are reconstructed and identified*/
	Int_t GetNumberOfIdentifiedEvents();
	/**The function loads the next event that passes the set cuts. It skips all events that do not fullfil the cuts. After event is loaded it values can be accessed by the Getter functions.
 	* @return true if a event with the specified cuts can be loaded, false otherwise. The cuts are defined by the parameters
 	* @param real - if true just events are loaded that are real events.
 	* @param iden - if true just events are loaded that are reconstructed and identified events.*/
 	
	Bool_t LoadEvent(Bool_t real,Bool_t iden);
	
	/**The function loads the event that passes the set cuts starting by the set event number. It skips all events that do not fullfil the cuts. After event is loaded it values can be accessed by the Getter functions.
 	* @return true if a event with the specified cuts can be loaded, false otherwise. The cuts are defined by the parameters
 	* @param no - number of the event in the tree where the loading starts.
 	* @param real - if true just events are loaded that are real events.
 	* @param iden - if true just events are loaded that are reconstructed and identified events.*/
	Bool_t LoadEvent(Int_t no,Bool_t real,Bool_t iden);
	/** Resets the counters in the class to load the tree entries in the file again*/
	void Reset();



	/**@return fPrimaryEn - The energy of the primary from the loaded real event*/
	Double_t GetRealPrimaryEnergy();
	/**@return fEnergy_p - The energy of the scattered photon from the loaded real event*/
	Double_t GetRealPhotonEnergy();
	/**@return fEnergy_e - The energy of the accelerated electron from the loaded real event*/
	Double_t GetRealElectronEnergy();
	/**@return *fPosition_p - The position of the first energy deposition of the scattered photon from the loaded real event*/
	TVector3 GetRealPhotonPosition();
	/**@return *fPosition_e - The position of the first energy deposition of the electron from the loaded real event*/
	TVector3 GetRealElectronPosition();
	/**@return *fPosition_s - The position of the primary photon origin from the loaded real event*/
	TVector3 GetRealSourcePosition();
	/**@return *fDirection_p - The direction of the scattered photon from the loaded real event*/
	TVector3 GetRealPhotonDirection();
	/**@return *fDirection_s - The direction of the primary photon leaving the source position from the loaded real event*/
	TVector3 GetRealPrimaryDirection();


	/**@return fReco_en_p->value - The energy of the scattered photon from the loaded reconstructed event*/
	Double_t GetRecoPhotonEnergy();
	/**@return fReco_en_p->value - The energy of the accelerated electron from the loaded reconstructed event*/
	Double_t GetRecoElectonEnergy();
	/**@return *fReco_pos_p->position - The position of the energy deposition of the scattered photon from the loaded reconstructed event*/
	TVector3 GetRecoPhotonPosition();
	/**@return *fReco_pos_e->position - The position of the energy deposition of the electron from the loaded reconstructed event*/
	TVector3 GetRecoElectronPosition();
	/**@return *fReco_dir_p->position - The direction of the scattered photon from the loaded reconstructed event*/
	TVector3 GetRecoPhotonDirection();

	/**@return fReco_en_p->uncertainty - The energy uncertainty of the scattered photon from the loaded reconstructed event*/
	Double_t GetRecoPhotonEnergyUn();
	/**@return fReco_en_p->uncertainty - The energy uncertainty of the accelerated electron from the loaded reconstructed event*/
	Double_t GetRecoElectonEnergyUn();
	/**@return *fReco_pos_p->uncertainty - The position uncertainty of the energy deposition of the scattered photon from the loaded reconstructed event*/
	TVector3 GetRecoPhotonPositionUn();
	/**@return *fReco_pos_e->uncertainty - The position uncertainty of the energy deposition of the electron from the loaded reconstructed event*/
	TVector3 GetRecoElectronPositionUn();
	/**@return *fReco_dir_p->uncertainty - The direction uncertainty of the scattered photon from the loaded reconstructed event*/
	TVector3 GetRecoPhotonDirectionUn();

  private:
	Bool_t verbose;//<! Indicates the verbose mode 

	Int_t fEventCo;//<! counter for the loaded events
	Int_t fEventCoReco;//<! counter for the loaded, identified events
	Int_t fEventCoReal;//<! counter for the loaded, real events
	
	TFile* fInputFile; //!< RootFile the input is save in 
	TTree* fInputTree_setup; //!< RootTree holding the CC setup information 
	TTree* fInputTree_real; //!< RootTree holding only the real event information 
	TTree* fInputTree_reco;//!< RootTree holding the reconstructed information 
	
	Int_t fEventNo;//!< eventNo of the events that are loaded from the file 

	TVector3* fPosition_scatterer;//!< fInputTree_setup BranchVariable - position of the scatterer
	TVector3* fPosition_absorber;//!< fInputTree_setup BranchVariable - position of the absorber 
	Double_t fScatterer_x; //!< fInputTree_setup BranchVariable - complete length of the scatterer in x
	Double_t fScatterer_y; //!< fInputTree_setup BranchVariable - complete length of the scatterer in y
	Double_t fScatterer_z; //!< fInputTree_setup BranchVariable - complete length of the scatterer in z
	Double_t fAbsorber_x; //!< fInputTree_setup BranchVariable - complete length of the absorber in x  
	Double_t fAbsorber_y; //!< fInputTree_setup BranchVariable - complete length of the absorber in y  
	Double_t fAbsorber_z; //!< fInputTree_setup BranchVariable - complete length of the absorber in z  

	
	Double_t fPrimaryEn; //!< fInputTree_real BranchVariable - real energy of the primary photon
	Double_t fEnergy_e; //!< fInputTree_real BranchVariable - real energy of electron
	Double_t fEnergy_p;//!< fInputTree_real BranchVariable - real energy of photon	
	TVector3* fPosition_e;	//!< fInputTree_real BranchVariable - real position of electron
	TVector3* fPosition_p;	//!< fInputTree_real BranchVariable - real position of photon
	TVector3* fPosition_s;	//!< fInputTree_real BranchVariable - source position
	TVector3* fDirection_p;//!< fInputTree_real BranchVariable - direction of photon after scattering
	TVector3* fDirection_s;//!< fInputTree_real BranchVariable - direction of primary leaving the source
                     
	Bool_t fIden;//!< states if the reco event was identified

	PhysicVar* fReco_en_e;//!< fInputTree_reco BranchVariable - reconstructed electron energy
	PhysicVar* fReco_en_p;//!< fInputTree_reco BranchVariable - reconstructed photon energy
	PhysicVec* fReco_pos_e;//!< fInputTree_reco BranchVariable - reconstructed electron position	
	PhysicVec* fReco_pos_p;//!< fInputTree_reco BranchVariable - reconstructed photon position	
	PhysicVec* fReco_dir_p;//!< fInputTree_reco BranchVariable - reconstructed photon direction
	std::vector<PhysicVec*>* fReco_cluster_pos;//!< fInputTree_reco BranchVariable - reconstructed cluster positions without classification
	std::vector<PhysicVar*>* fReco_cluster_en;//!< fInputTree_reco BranchVariable - reconstructed cluister energies without classification

};


#endif //G4INPUT_H
