#include "G4Input.hh"

G4Input::G4Input(TString InputFile, Bool_t verb)
    : verbose(verb), fEventCo(0), fEventCoReal(0), fEventCoReco(0), fEventNo(0), fPrimaryEn(0),
      fEnergy_e(0), fEnergy_p(0), fIden(false), fScatterer_x(0), fScatterer_y(0), fScatterer_z(0),
      fAbsorber_x(0), fAbsorber_y(0), fAbsorber_z(0)
{
    fPosition_scatterer = new TVector3(0, 0, 0);
    fPosition_absorber = new TVector3(0, 0, 0);

    fPosition_e = new TVector3(0, 0, 0);
    fPosition_p = new TVector3(0, 0, 0);
    fPosition_s = new TVector3(0, 0, 0);
    fDirection_p = new TVector3(0, 0, 0);
    fDirection_s = new TVector3(0, 0, 0);

    fReco_en_e = new PhysicVar(0, 0);
    fReco_en_p = new PhysicVar(0, 0);
    fReco_pos_e = new PhysicVec(TVector3(0, 0, 0), TVector3(0, 0, 0));
    fReco_pos_p = new PhysicVec(TVector3(0, 0, 0), TVector3(0, 0, 0));
    fReco_dir_p = new PhysicVec(TVector3(0, 0, 0), TVector3(0, 0, 0));

    fReco_cluster_pos = 0;
    fReco_cluster_en = 0;

    if (InputFile != "")
    {
        fInputFile = new TFile(InputFile, "READ");

        if (fInputFile->IsZombie() != true)
        {
            fInputTree_setup = (TTree*)fInputFile->Get("G4SimulationData_Setup");

            fInputTree_setup->SetBranchAddress("ScattererThickness_x", &fScatterer_x);
            fInputTree_setup->SetBranchAddress("ScattererThickness_y", &fScatterer_y);
            fInputTree_setup->SetBranchAddress("ScattererThickness_z", &fScatterer_z);
            fInputTree_setup->SetBranchAddress("AbsorberThickness_x", &fAbsorber_x);
            fInputTree_setup->SetBranchAddress("AbsorberThickness_y", &fAbsorber_y);
            fInputTree_setup->SetBranchAddress("AbsorberThickness_z", &fAbsorber_z);
            fInputTree_setup->SetBranchAddress("ScattererPosition", &fPosition_scatterer);
            fInputTree_setup->SetBranchAddress("AbsorberPosition", &fPosition_absorber);
            fInputTree_setup->GetEntry();

            fInputTree_real = (TTree*)fInputFile->Get("G4SimulationData_Real");

            fInputTree_real->SetBranchAddress("EventNumber", &fEventNo);
            fInputTree_real->SetBranchAddress("Energy_Primary", &fPrimaryEn);
            fInputTree_real->SetBranchAddress("RealEnergy_e", &fEnergy_e);
            fInputTree_real->SetBranchAddress("RealEnergy_p", &fEnergy_p);
            fInputTree_real->SetBranchAddress("RealPosition_source", &fPosition_s);
            fInputTree_real->SetBranchAddress("RealDirection_source", &fDirection_s);
            fInputTree_real->SetBranchAddress("RealPosition_e", &fPosition_e);
            fInputTree_real->SetBranchAddress("RealPosition_p", &fPosition_p);
            fInputTree_real->SetBranchAddress("RealDirection_scatter", &fDirection_p);

            fInputTree_reco = (TTree*)fInputFile->Get("G4SimulationData_Reconstruction");

            fInputTree_reco->SetBranchAddress("EventNumber", &fEventNo);
            fInputTree_reco->SetBranchAddress("Identified", &fIden);
            fInputTree_reco->SetBranchAddress("RecoEnergy_e", &fReco_en_e);
            fInputTree_reco->SetBranchAddress("RecoEnergy_p", &fReco_en_p);
            fInputTree_reco->SetBranchAddress("RecoPosition_e", &fReco_pos_e);
            fInputTree_reco->SetBranchAddress("RecoPosition_p", &fReco_pos_p);
            fInputTree_reco->SetBranchAddress("RecoDirection_scatter", &fReco_dir_p);
            fInputTree_reco->SetBranchAddress("RecoClusterPositions", &fReco_cluster_pos);
            fInputTree_reco->SetBranchAddress("RecoClusterEnergies", &fReco_cluster_en);
        }
        else
        {
            std::cout << "The inputfile with the name " << InputFile << " is zombie! Abort!"
                      << std::endl;
            std::abort();
        }
    }
    else
    {
        std::cout << "The name of the inputfile for the G4Input is not defined! Abort!"
                  << std::endl;
        std::abort();
    }
}

G4Input::~G4Input()
{
    delete fPosition_scatterer;
    delete fPosition_absorber;
    delete fPosition_e;
    delete fPosition_p;
    delete fPosition_s;
    delete fDirection_p;
    delete fDirection_s;
    delete fReco_en_e;
    delete fReco_en_p;
    delete fReco_pos_e;
    delete fReco_pos_p;
    delete fReco_dir_p;
    delete fReco_cluster_pos;
    delete fReco_cluster_en;
}

Bool_t G4Input::LoadEvent(Bool_t real = false, Bool_t iden = false)
{
    if (real == false)
    {
        if (iden == false) { return LoadEvent(fEventCo, real, iden); }
        else
        {
            return LoadEvent(fEventCoReco, real, iden);
        }
    }
    else
    {
        return LoadEvent(fEventCoReal, real, iden);
    }
}

Bool_t G4Input::LoadEvent(Int_t no, Bool_t real = false, Bool_t iden = false)
{
    Bool_t error = true;
    if (real == false)
    {
        if (iden == false)
        {
            fInputTree_real->GetEntry(no);
            fInputTree_reco->GetEntry(no);
            error = false;
            fEventCo = no + 1;
        }
        else
        {
            for (int i = no; i < fInputTree_reco->GetEntries(); i++)
            {
                fInputTree_reco->GetEntry(i);
                if (fIden != false)
                {
                    if (verbose)
                    {
                        std::cout << "Identified is requested. The loaded event has the "
                                     "event number "
                                  << fEventNo << "." << std::endl;
                    }
                    fInputTree_real->GetEntry(i);
                    fEventCoReco = i + 1;
                    error = false;
                    break;
                }
            }
        }
    }
    else
    {
        for (int i = no; i < fInputTree_real->GetEntries(); i++)
        {
            fInputTree_real->GetEntry(i);
            if (fEnergy_e != 0)
            {
                if (iden == false)
                {
                    if (verbose)
                    {
                        std::cout << "Real event is requested. The loaded event has the "
                                     "event number "
                                  << fEventNo << "." << std::endl;
                    }
                    fInputTree_reco->GetEntry(i);
                    fEventCoReal = i + 1;
                    error = false;
                    break;
                }
                else
                {
                    fInputTree_reco->GetEntry(i);
                    if (fIden != false)
                    {
                        if (verbose)
                        {
                            std::cout << "Real event and identified is requested. The loaded "
                                         "event has the event number "
                                      << fEventNo << "." << std::endl;
                        }
                        fEventCoReal = i + 1;
                        error = false;
                        break;
                    }
                }
            }
        }
    }
    if (error)
    {
        if (no == 0) std::cerr << "No event with the requested specification found!" << std::endl;
        return false;
    }
    else
        return true;
}

void G4Input::Reset()
{
    if (verbose) { std::cout << "The G4Fileloader is reseted!" << std::endl; }
    fEventCo = 0;
    fEventCoReal = 0;
    fEventCoReco = 0;
}

Int_t G4Input::GetNumberOfRealEvents() { return fInputTree_real->GetEntries("RealEnergy_e!=0"); }

Int_t G4Input::GetNumberOfEvents() { return fInputTree_reco->GetEntries(); }

Int_t G4Input::GetNumberOfIdentifiedEvents()
{
    return fInputTree_real->GetEntries("Identified==true");
}

// Getter for setup

Double_t G4Input::GetScattererXLength() { return fScatterer_x; }

Double_t G4Input::GetScattererYLength() { return fScatterer_y; }

Double_t G4Input::GetScattererZLength() { return fScatterer_z; }

Double_t G4Input::GetAbsorberXLength() { return fAbsorber_x; }

Double_t G4Input::GetAbsorberYLength() { return fAbsorber_y; }

Double_t G4Input::GetAbsorberZLength() { return fAbsorber_z; }

TVector3 G4Input::GetScattererPosition() { return *fPosition_scatterer; }

TVector3 G4Input::GetAbsorberPosition() { return *fPosition_absorber; }

// Getter for the real event

Double_t G4Input::GetRealPrimaryEnergy() { return fPrimaryEn; }

Double_t G4Input::GetRealPhotonEnergy() { return fEnergy_p; }

Double_t G4Input::GetRealElectronEnergy() { return fEnergy_e; }

TVector3 G4Input::GetRealPhotonPosition() { return *fPosition_p; }

TVector3 G4Input::GetRealSourcePosition() { return *fPosition_s; }

TVector3 G4Input::GetRealPhotonDirection() { return *fDirection_p; }

TVector3 G4Input::GetRealPrimaryDirection() { return *fDirection_s; }

TVector3 G4Input::GetRealElectronPosition() { return *fPosition_e; }

// Getter reco event

Double_t G4Input::GetRecoPhotonEnergy() { return fReco_en_p->value; }

Double_t G4Input::GetRecoElectonEnergy() { return fReco_en_e->value; }

TVector3 G4Input::GetRecoPhotonPosition() { return fReco_pos_p->position; }

TVector3 G4Input::GetRecoElectronPosition() { return fReco_pos_e->position; }

TVector3 G4Input::GetRecoPhotonDirection() { return fReco_dir_p->position; }

Double_t G4Input::GetRecoPhotonEnergyUn() { return fReco_en_p->uncertainty; }

Double_t G4Input::GetRecoElectonEnergyUn() { return fReco_en_e->uncertainty; }

TVector3 G4Input::GetRecoPhotonPositionUn() { return fReco_pos_p->uncertainty; }

TVector3 G4Input::GetRecoElectronPositionUn() { return fReco_pos_e->uncertainty; }

TVector3 G4Input::GetRecoPhotonDirectionUn() { return fReco_dir_p->uncertainty; }
