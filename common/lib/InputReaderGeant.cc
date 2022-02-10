#include "InputReaderGeant.hh"
#include <vector>
ClassImp(InputReaderGeant);

//------------------------------------------------------------------
InputReaderGeant::InputReaderGeant() : InputReader() {
  Clear();
  cout << "##### Warning in InputReaderGeant constructor!" << endl;
  cout << "You are usinf default constructor." << endl;
}
//------------------------------------------------------------------
InputReaderGeant::InputReaderGeant(TString path):
InputReader(path),
fLoadReal(false),
fCorrectOnly(0),
fPrimaryEnergy(0),
fEnergyLoss(0), 
fEnergyScattered(0) 
{
  if (!AccessTree("Events","Setup")) {
    throw "##### Exception in InputReaderGeant constructor!";
  }
  fPositionScat = new TVector3();
  fPositionAbs = new TVector3();
  fDirectionScat = new TVector3();
  
  fPositionSource = new TVector3();
  fDirectionSource = new TVector3();


}
//------------------------------------------------------------------
InputReaderGeant::~InputReaderGeant() {
  if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
bool InputReaderGeant::AccessTree(TString name, TString namesetup) {

  fTree = (TTree*)fFile->Get(name);
  fTreeSetup = (TTree*)fFile->Get(namesetup);
  //fTree2 = (TTree*)fFile->Get(name2);
  if (fTree == NULL) {
    cout << "##### Error in InputReaderGeant::AccessTree()!" << endl;
    cout << "Could not access the tree with requested name" << name.Data() << endl;
    return false;
  }
  if (fTreeSetup == NULL) {
    cout << "##### Error in InputReaderGeant::AccessTree()!" << endl;
    cout << "Could not access the tree with requested name" << namesetup.Data() << endl;
    return false;
  }
  
  fRecoEnergy_e = new PhysicVar();
  fRecoEnergy_p = new PhysicVar();
  fRecoPosition_e = new PhysicVec();
  fRecoPosition_p = new PhysicVec();
  fRecoDirection_scatter = new PhysicVec();
  fRecoClusterPositions = new vector<PhysicVec>;
  fRecoClusterEnergies = new vector<PhysicVar>;
  
  fRealPosition_source = new TVector3();
  fRealDirection_source = new TVector3();
  fRealComptonPosition = new TVector3();
  fRealDirection_scatter = new TVector3();
  
  fRealPosition_e = new vector<TVector3>;
  fRealInteractions_e = new vector<int>;
  fRealPosition_p = new vector<TVector3>;
  fRealInteractions_p = new vector<int>;
  
///////////////////////////////////////////////  
  
  fScattererPosition = new TVector3();
  fAbsorberPosition = new TVector3();
  
  fScattererDim = new TVector3();
  fAbsorberDim = new TVector3();
 
  double scattererThickness_x=0; 
  double scattererThickness_y=0; 
  double scattererThickness_z=0; 
  double absorberThickness_x=0;
  double absorberThickness_y=0;
  double absorberThickness_z=0;
 
  fTreeSetup->SetBranchAddress("ScattererThickness_x", &scattererThickness_x);
  fTreeSetup->SetBranchAddress("ScattererThickness_y", &scattererThickness_y);
  fTreeSetup->SetBranchAddress("ScattererThickness_z", &scattererThickness_z);
  fTreeSetup->SetBranchAddress("AbsorberThickness_x", &absorberThickness_x);
  fTreeSetup->SetBranchAddress("AbsorberThickness_y", &absorberThickness_y);
  fTreeSetup->SetBranchAddress("AbsorberThickness_z", &absorberThickness_z);
  fTreeSetup->SetBranchAddress("ScattererPosition", &fScattererPosition);
  fTreeSetup->SetBranchAddress("AbsorberPosition", &fAbsorberPosition);
  fTreeSetup->SetBranchAddress("NumberOfSimulatedEvents",
                           &fNumberOfSimulatedEvents);
  
  fTree->SetBranchAddress("EventNumber", &fEventNumber);
  fTree->SetBranchAddress("Energy_Primary", &fEnergy_Primary);
  fTree->SetBranchAddress("SimulatedEventType", &fSimulatedEventType);
  fTree->SetBranchAddress("RealEnergy_e", &fRealEnergy_e);
  fTree->SetBranchAddress("RealEnergy_p", &fRealEnergy_p);
  fTree->SetBranchAddress("RealPosition_source", &fRealPosition_source);
  fTree->SetBranchAddress("RealDirection_source", &fRealDirection_source);
  fTree->SetBranchAddress("RealPosition_e", &fRealPosition_e);
  fTree->SetBranchAddress("RealInteractions_e", &fRealInteractions_e);
  fTree->SetBranchAddress("RealComptonPosition", &fRealComptonPosition);
  fTree->SetBranchAddress("RealPosition_p", &fRealPosition_p);
  fTree->SetBranchAddress("RealInteractions_p", &fRealInteractions_p);
  fTree->SetBranchAddress("RealDirection_scatter", &fRealDirection_scatter);
  fTree->SetBranchAddress("Identified", &fIdentified);
  fTree->SetBranchAddress("PurCrossed", &fPurCrossed);
  fTree->SetBranchAddress("RecoEnergy_e", &fRecoEnergy_e);
  fTree->SetBranchAddress("RecoEnergy_p", &fRecoEnergy_p);
  fTree->SetBranchAddress("RecoPosition_e", &fRecoPosition_e);
  fTree->SetBranchAddress("RecoPosition_p", &fRecoPosition_p);
  fTree->SetBranchAddress("RecoDirection_scatter", &fRecoDirection_scatter);
  fTree->SetBranchAddress("RecoClusterPositions", &fRecoClusterPositions);
  fTree->SetBranchAddress("RecoClusterEnergies", &fRecoClusterEnergies);

  cout << "\n\nIn InputReaderGeant::AccessTree()." << endl;
  cout << fTree->GetName() << " tree accessed.\n" << endl;

  cout << "\n\nIn InputReaderGeant::AccessTree()." << endl;
  cout << fTreeSetup->GetName() << " tree accessed.\n" << endl;
  
  fTreeSetup->GetEntry(0);
  fScattererDim->SetXYZ(scattererThickness_x,scattererThickness_y,scattererThickness_z);
  fAbsorberDim->SetXYZ(absorberThickness_x,absorberThickness_y,absorberThickness_z);
  return true;
}
void InputReaderGeant::SelectEvents(){  
  fSelectedEvents.clear(); 
  for(int i=0;i<fTree->GetEntries();i++) {
  	fTree->GetEntry(i);
 	if(SelectSingleEvent())fSelectedEvents.push_back(i);
  }
}	
bool InputReaderGeant::SelectSingleEvent(){
	bool returner=true;
	if(fLoadReal){
  	      // For the detector also only events where interactions in both detector moduls occure are triggered, for better comparison also only these real events are used. 
  	      bool physinteractions[2]={false,false};
  	      for(auto k=0;k<fRecoClusterPositions->size();k++){
  	              if(fRecoClusterPositions->at(k).position.x()<(fScattererPosition->X()+fScattererDim->X()/2))physinteractions[0]=true;
  	              else if(fRecoClusterPositions->at(k).position.x()>(fAbsorberPosition->X()-fAbsorberDim->X()/2))physinteractions[1]=true;
  	      }
  	      if(!physinteractions[0] || !physinteractions[1] ) return false;
  	      // Only Compton Events with at least one additional interaction of the photon( first one is the compton effect (fRealInteraction.at(0) is the compton effect itself) and at least one interaction of the electron can be reconstructed.
  	      if(fSimulatedEventType!=2) return false;
  	      if(!(fRealEnergy_e>0)) return false;
  	      if(fRealInteractions_p->size() < 2 || fRealInteractions_e->size() < 1) return false; 
  	      // Also a selection with one of the particles interacting in both modules (besides the first compton effect) are difficult to select so these ones are not chosen as well
  	      //indicates if photon/electron inetacted in scatterer/absorber
  	      //photonscatter/electronscatterer/photonabsorber/electronscatterer
              bool interactions[4]={false,false,false,false};
	      bool photosec=false;
	      bool electronsec=false;
  	      for(int m = 1; m < fRealInteractions_p->size(); m++) {
		if(fRealInteractions_p->at(m)<10 && fRealInteractions_p->at(m)>0)photosec=true;
  	      	if(fRealPosition_p->at(m).X()<(fScattererPosition->X()+fScattererDim->X()/2))interactions[0]=true;
  	      	else if(fRealPosition_p->at(m).X()>(fAbsorberPosition->X()-fAbsorberDim->X()/2))interactions[2]=true;
  	      }
  	      for(int m = 0; m < fRealInteractions_e->size(); m++) {
		if(fRealInteractions_e->at(m)<20)electronsec=true;
  	      	if(fRealPosition_e->at(m).X()<(fScattererPosition->X()+fScattererDim->X()))interactions[1]=true;
  	      	else if(fRealPosition_e->at(m).X()>(fAbsorberPosition->X()-fAbsorberDim->X()))interactions[3]=true;
  	      }
  	      if(!photosec || !electronsec) return false;
              if(interactions[0]&&interactions[1])return false;
  	      if(interactions[2]&&interactions[3])return false;
  	}
  	else{
  	      if(fIdentified==0) return false;
  	      if(fCorrectOnly){
  	        if(!(fRealEnergy_e>0)) returner=false;
  	        if(fSimulatedEventType!=2) returner=false;
        	if(fRealInteractions_p->size() < 2 || fRealInteractions_e->size() < 1)returner=false; 
		bool interactions[4]={false,false,false,false};
	      	bool photosec=false;
	      	bool electronsec=false;
		int e_pos_index=-1;
		int p_pos_index=-1;
        	for(int m = 1; m < fRealInteractions_p->size(); m++) {
			if(fRealInteractions_p->at(m)<10 && fRealInteractions_p->at(m)>0){
				if(photosec==false)p_pos_index=m;	
				photosec=true;
			}
        		if(fRealPosition_p->at(m).X()<(fScattererPosition->X()+fScattererDim->X()))interactions[0]=true;
        		else if(fRealPosition_p->at(m).X()>(fAbsorberPosition->X()-fAbsorberDim->X()))interactions[2]=true;
		}
        	for(int m = 0; m < fRealInteractions_e->size(); m++) {
			if(fRealInteractions_e->at(m)<20){
				if(electronsec==false)e_pos_index=m;	
				electronsec=true;
			}
        		if(fRealPosition_e->at(m).X()<(fScattererPosition->X()+fScattererDim->X()))interactions[1]=true;
        		else if(fRealPosition_e->at(m).X()>(fAbsorberPosition->X()-fAbsorberDim->X()))interactions[3]=true;
        	}
  	      	if(!photosec || !electronsec) returner=false;
		if(interactions[0]==interactions[1])returner=false;
		if(interactions[2]==interactions[3])returner=false;
		//checking the asked events
		if(fCorrectOnly==3 && returner==false){
			return true;
		}
		else if (fCorrectOnly==3)returner=false;
		if(returner){
			if(fCorrectOnly==1){
				return true;
			}
			if((fCorrectOnly==2 || fCorrectOnly==4)){
				//Compare reco and real
				TVector3 uncvec(2.6,30,2.6);
				double uncen=0.12;
				TVector3 diffvec;
				double diffen=0;
				double en=0;
				for(int u=0;u<2;u++){
					if(u==0){
						diffvec= fRealPosition_e->at(e_pos_index)-fRecoPosition_e->position;
						diffen=fRealEnergy_e-fRecoEnergy_e->value; 
						en=fRecoEnergy_e->value;
					}
					else{
						diffvec= fRealPosition_p->at(p_pos_index)-fRecoPosition_p->position;
						diffen=fRealEnergy_p-fRecoEnergy_p->value; 
						en=fRecoEnergy_p->value;
					}
					if(fabs(diffvec.x())>(uncvec.x()))returner=false;
					if(fabs(diffvec.y())>(uncvec.y()))returner=false;
					if(fabs(diffvec.z())>(uncvec.z()))returner=false;
					if((fabs(diffen)/en)>uncen) returner=false;
				}
				if(fCorrectOnly==2 && returner==true) return true;
				if(fCorrectOnly==4 && returner==false) return true;
				else if(fCorrectOnly==4) returner=false;
			}
		}
	      }
	}
	return returner;
}
//------------------------------------------------------------------
bool InputReaderGeant::LoadEvent(int i) {

  int imax = fTree->GetEntries();
  if (i > imax) {
    cout << "##### Error in InputReaderGeant::LoadEvent() in Event tree!" << endl;
    cout << "Requested event number larger than number of events in the tree!"
         << endl;
    return false;
  }
  fTree->GetEntry(i);
  //Freakycombo
  if(SelectSingleEvent()==false)return false;
  if(fLoadReal){
  	for(int m = 1; m < fRealInteractions_p->size(); m++) {
		if(fRealInteractions_p->at(m)<10 &&fRealInteractions_p->at(m)>0){
       			fPositionAbs->SetXYZ(fRealPosition_p->at(m).X(), fRealPosition_p->at(m).Y(), fRealPosition_p->at(m).Z());
			break;
		}
	}
  	for(int m = 1; m < fRealInteractions_e->size(); m++) {
		if(fRealInteractions_e->at(m)<20){
       			fPositionScat->SetXYZ(fRealPosition_e->at(m).X(), fRealPosition_e->at(m).Y(), fRealPosition_e->at(m).Z());
			break;
		}
	}
  	fDirectionScat->SetXYZ(fRealDirection_scatter->X(),fRealDirection_scatter->Y(),fRealDirection_scatter->Z());
  	
	fPrimaryEnergy = fRealEnergy_e + fRealEnergy_p;
	fEnergyLoss=fRealEnergy_e; 
	fEnergyScattered=fRealEnergy_p; 
  }
  else{
	fPositionScat->SetXYZ(fRecoPosition_e->position.X(),fRecoPosition_e->position.Y(),fRecoPosition_e->position.Z());
  	fPositionAbs->SetXYZ(fRecoPosition_p->position.X(),fRecoPosition_p->position.Y(),fRecoPosition_p->position.Z());
	fDirectionScat->SetXYZ(fRecoDirection_scatter->position.X(),fRecoDirection_scatter->position.Y(),fRecoDirection_scatter->position.Z());
  	fPrimaryEnergy = fRecoEnergy_e->value + fRecoEnergy_p->value;
	fEnergyLoss=fRecoEnergy_e->value; 
	fEnergyScattered=fRecoEnergy_p->value; 
  }

  return true;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionPrimary(void) {
  if(fLoadReal){
  	cout << "##### Warning in InputReaderGeant::GetPositionPrimary()!" << endl;
  	cout << "\t Position of gamma source is unknown!" << endl;
  	return NULL;
  }
  else{
  	return fRealPosition_source;
  }
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetGammaDirPrimary(void) {
  if(fLoadReal){
  	cout << "##### Warning in InputReaderGeant::GetGammaDirPrimary()!" << endl;
  	cout << "\t Direction of primary gamma is unknown!" << endl;
  	return NULL;
  }
  else{
  	return fDirectionSource;
  }
}
void InputReaderGeant::Clear(void) {
  fEventNumber = -1;
  fIdentified = -1000;
  fPurCrossed = false;
  fSimulatedEventType= -1000;
  fEnergy_Primary = -1000;
  fRealEnergy_e = -1000;
  fRealEnergy_p = -1000;
  fRealPosition_e = NULL;
  fRealComptonPosition = NULL;
  fRealPosition_p = NULL;
  fPositionScat = NULL;
  fPositionAbs = NULL;
  fDirectionScat = NULL;
  fDirectionSource = NULL;
  fPositionSource = NULL;
  fTreeSetup = NULL;
  fTree = NULL;
  fFile = NULL;
  fRealDirection_scatter = NULL;
  fRealDirection_source = NULL;
  fRealPosition_source = NULL;
  
  fRecoEnergy_e = NULL;
  fRecoEnergy_p = NULL;
  fRecoPosition_e = NULL;
  fRecoPosition_p = NULL;
  
  if (!fRecoClusterPositions->empty()) fRecoClusterPositions->clear();
  if (!fRecoClusterEnergies->empty()) fRecoClusterEnergies->clear();
  
  fNumberOfSimulatedEvents = -1;

  fScattererPosition = NULL;
  fAbsorberPosition = NULL;
  
  fScattererDim= NULL;
  fAbsorberDim= NULL;

  fLoadReal=false;
  fPrimaryEnergy=-1000;
  fEnergyLoss=-1000;
  fEnergyScattered=-1000;

  fPositionScat = NULL;
  fPositionAbs = NULL;
  fDirectionScat = NULL;

  return;
}
//------------------------------------------------------------------
