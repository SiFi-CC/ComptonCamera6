#include <TString.h>
#include <TApplication.h>
#include <TH1D.h>
#include <TCanvas.h>

#include "G4Input.hh"

int main(int argc,char** argv){
	TString name ="test.root";

	TApplication* TApp= new TApplication("TheApp",0,0);

	G4Input* test = new G4Input(name,true);

	//Example how to get the setup stuff
	
	std::cout << "The scatterer is positioned at " << test->GetScattererPosition().x() << ", " << test->GetScattererPosition().y() << ", " << test->GetScattererPosition().z() << " and has the Dimensions: x: " << test->GetScattererXLength() << " y: " << test->GetScattererYLength() << " z: " <<test->GetScattererZLength() <<std::endl; 
	std::cout << "The absorber is positioned at " << test->GetAbsorberPosition().x() << ", " << test->GetAbsorberPosition().y() << ", " << test->GetAbsorberPosition().z() << " and has the Dimensions: x: " << test->GetAbsorberXLength() << " y: " << test->GetAbsorberYLength() << " z: " <<test->GetAbsorberZLength() <<std::endl; 

	//Example how to get the event variables

	TH1D* Energy_P_real_all = new TH1D("EPRA","EPRA",700,0,7); //Contains the Energy of the real Photons without any cuts
	TH1D* Energy_P_real_cut = new TH1D("EPRC","EPRC",700,0,7);//Contains the Energy of the real Photons cut on events where actually a compton scattering occured
	TH1D* Energy_P_reco_all = new TH1D("EPReA","EPReA",700,0,7);//Contains the Energy of the reconstructed Photons without any cuts
	TH1D* Energy_P_reco_cut = new TH1D("EPReC","EPReC",700,0,7);//Contains the Energy of the reconstructed Photons cut on events that are identified
	TH1D* Energy_P_real_reco= new TH1D("EPRR","EPRR",700,0,7);//Contains the Energy of the reconstructed Photons cut on events that are real events and identified
	//to catch in generall all in the input file stored events (means reconstructed at all)
	for(int i=0;i<test->GetNumberOfEvents();i++){
		if(test->LoadEvent(false,false)==false) break;
		Energy_P_real_all->Fill(test->GetRealPhotonEnergy());	
		Energy_P_reco_all->Fill(test->GetRecoPhotonEnergy());	
	}
	//to catch all in the input file stored events that contain a real compton scattering proccess in the detector
	for(int i=0;i<test->GetNumberOfEvents();i++){
		if(test->LoadEvent(true,false)==false) break;
		Energy_P_real_cut->Fill(test->GetRealPhotonEnergy());	
		std::cout << "The real photon position is at " <<  test->GetRealPhotonPosition().x() << ", " <<  test->GetRealPhotonPosition().y() << ", " <<  test->GetRealPhotonPosition().z() << std::endl; 
	}
	//to catch all in the input file stored events that are reconstructed and identified
	for(int i=0;i<test->GetNumberOfEvents();i++){
		if(test->LoadEvent(false,true)==false)break;
		Energy_P_reco_cut->Fill(test->GetRecoPhotonEnergy());	
		std::cout << "The reco photon position is at " <<  test->GetRecoPhotonPosition().x() << ", " <<  test->GetRecoPhotonPosition().y() << ", " <<  test->GetRecoPhotonPosition().z() << std::endl; 
	}
	//to catch all in the input file stored events that are reconstructed and identified and also a compton scattering occured in the detector
	for(int i=0;i<test->GetNumberOfEvents();i++){
		if(test->LoadEvent(true,true)==false) break;
		Energy_P_real_reco->Fill(test->GetRecoPhotonEnergy());	
	}

	TCanvas* tempCan= new TCanvas("temp","temp",1600,1600);
	tempCan->Divide(2,3);
	tempCan->cd(1);
	Energy_P_real_all->Draw("");
	tempCan->cd(2);
	Energy_P_real_cut->Draw("");
	tempCan->cd(3);
	Energy_P_reco_all->Draw("");
	tempCan->cd(4);
	Energy_P_reco_cut->Draw("");
	tempCan->cd(5);
	Energy_P_real_reco->Draw("");
	//To loop again over the input file 

	test->Reset();
	
	//to catch all in the input file stored events that are reconstructed and identified
	for(int i=0;i<test->GetNumberOfEvents();i++){
		if(test->LoadEvent(false,true)==false)break;
		std::cout << "The reco electron position is at " <<  test->GetRecoElectronPosition().x() << ", " <<  test->GetRecoElectronPosition().y() << ", " <<  test->GetRecoElectronPosition().z() << std::endl; 
	}


	TApp->Run();

}
