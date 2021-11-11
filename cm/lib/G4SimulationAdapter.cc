#include "G4SimulationAdapter.hh"
#include <TFile.h>
#include <TParameter.h>
#include <TSystem.h>

#include "CmdLineConfig.hh"

G4SimulationAdapter::G4SimulationAdapter(TString filename)
{
    fSelected = new TFile(filename, "READ");
}

G4SimulationAdapter::~G4SimulationAdapter()
{
    for (auto file : fFiles)
    {
        delete file;
    }
}

std::vector<TFile*> G4SimulationAdapter::Filter(std::function<bool(TFile*)> filter)
{
    std::vector<TFile*> filtered;
    std::copy_if(fSimulations.begin(), fSimulations.end(), std::back_inserter(filtered), filter);
    return filtered;
}

/**
 *  Potentially file could contain data for multiple simulations, e.g. different
 * energies, distance to mask, distance to detector, material etc. This function
 * will pick whathever is first on the list and find all matching to it.
 */
CameraGeometry G4SimulationAdapter::GetFirstReconstructData()
{
    CameraGeometry geometry;
    ParseSelected(&geometry);
    fSelected->cd();
    geometry.fMatrixHCam.Read("matrixH");

    return geometry;
}

bool G4SimulationAdapter::IsSimulationGeometryEqual(TList* sim1, TList* sim2)
{
    /* very inefficient implementation */
    for (auto entryFromFirst : *sim1)
    {
        if (TString(entryFromFirst->ClassName()) != "TParameter<double>") { continue; }
        auto e1 = static_cast<TParameter<double>*>(entryFromFirst);
        if (fIgnoredKeys.find(std::string(e1->GetName())) != fIgnoredKeys.end()) { continue; }

        for (auto entry : *sim2)
        {
            if (TString(entry->ClassName()) != "TParameter<double>") { continue; }
            auto e2 = static_cast<TParameter<double>*>(entry);
            if (TString(e1->GetName()) == TString(e2->GetName()) && e1->GetVal() != e2->GetVal())
            {
                log->debug("skipping because key={} has different values ({}, {})", e1->GetName(),
                           e1->GetVal(), e2->GetVal());
                return false;
            }
        }
    }
    return true;
}

void G4SimulationAdapter::VerifyForReconstruct(TFile* simulationFile)
{
    auto firstUserInfo = static_cast<TTree*>(fSelected->Get("metadata"))->GetUserInfo();
    auto detector = static_cast<TTree*>(simulationFile->Get("metadata"))->GetUserInfo();

    if (!IsSimulationGeometryEqual(detector, firstUserInfo))
    {
        log->error("Geometry of detector simulation does not match geometry of provided "
                   "grid simulation data");
        // throw "input data mismatch";
    }
}

void G4SimulationAdapter::ParseSelected(CameraGeometry* camera)
{
    auto userInfo = static_cast<TTree*>(fSelected->Get("metadata"))->GetUserInfo();
    std::unordered_map<std::string, double> map;
    for (auto entry : *userInfo)
    {
        if (TString(entry->ClassName()) != "TParameter<double>") { continue; }
        auto e = static_cast<TParameter<double>*>(entry);
        map[std::string(e->GetName())] = e->GetVal();
    }

    camera->detector.xRange = {map["detectorMinX"], map["detectorMaxX"]};
    camera->detector.yRange = {map["detectorMinY"], map["detectorMaxY"]};
    camera->detector.zRange = {map["detectorMinZ"], map["detectorMaxZ"]};
    camera->detector.binX = map["detectorBinX"];
    camera->detector.binY = map["detectorBinY"];
    camera->detector.binZ = 1;

    camera->mask.xRange = {map["maskMinX"], map["maskMaxX"]};
    camera->mask.yRange = {map["maskMinY"], map["maskMaxY"]};
    camera->mask.zRange = {map["maskMinZ"], map["maskMaxZ"]};
    camera->mask.binX = map["maskBinX"];
    camera->mask.binY = map["maskBinY"];
    camera->mask.binZ = 1;

    camera->source.xRange = {map["sourceMinX"], map["sourceMaxX"]};
    camera->source.yRange = {map["sourceMinY"], map["sourceMaxY"]};
    camera->source.zRange = {-0.1, 0.1};
    camera->source.binX = map["sourceBinX"];
    camera->source.binY = map["sourceBinY"];
    camera->source.binZ = 1;
}
