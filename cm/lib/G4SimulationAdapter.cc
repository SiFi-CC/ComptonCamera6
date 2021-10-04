#include "G4SimulationAdapter.hh"

#include <TFile.h>
#include <TKey.h>
#include <TParameter.h>
#include <TSystem.h>

G4SimulationAdapter::G4SimulationAdapter(TString filename)
{
    fFiles.push_back(new TFile(filename, "READ"));
    log->info("{}", TString::Format("%s.%d", filename.Data(), 1).Data());
    for (int i = 1; !gSystem->AccessPathName(TString::Format("%s.%d", filename.Data(), i)); i++)
    {
        fFiles.push_back(new TFile(TString::Format("%s.%d", filename.Data(), i), "READ"));
    }
    for (auto file : fFiles)
    {
        if (!file->IsOpen())
        {
            log->error("Unable to open file {}", file->GetName());
            throw "unable to open file";
        }
    }
    ReadMetadata();
}

G4SimulationAdapter::~G4SimulationAdapter()
{
    for (auto file : fFiles)
    {
        delete file;
    }
}

void G4SimulationAdapter::ReadMetadata()
{
    for (auto file : fFiles)
    {
        if (file->IsFolder())
        {
            log->info("Found {} data files", file->GetListOfKeys()->GetSize());
        }
        auto fileContent = file->GetListOfKeys();
        for (auto element : *fileContent)
        {
            auto el = static_cast<TKey*>(element);
            fSimulations.push_back(static_cast<TFile*>(file->Get(el->GetName())));
        }
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
    fSelected = *fSimulations.begin();
    auto firstUserInfo = static_cast<TTree*>(fSelected->Get("metadata"))->GetUserInfo();

    log->info("Picked reconstruction data for:");
    for (auto entryFromFirst : *firstUserInfo)
    {
        if (TString(entryFromFirst->ClassName()) != "TParameter<double>") { continue; }
        auto e = static_cast<TParameter<double>*>(entryFromFirst);
        if (fIgnoredKeys.find(std::string(e->GetName())) != fIgnoredKeys.end()) { continue; }
        log->info("key={}, value={}", e->GetName(), e->GetVal());
    }

    CameraGeometry geometry;
    ParseSelected(&geometry);
    geometry.recoData = Filter(
        [firstUserInfo, this](TFile* file) -> bool
        {
            auto userInfo = static_cast<TTree*>(file->Get("metadata"))->GetUserInfo();
            return IsSimulationGeometryEqual(userInfo, firstUserInfo);
        });
    ParsePointSources(&geometry);
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
        throw "input data mismatch";
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
}

void G4SimulationAdapter::ParsePointSources(CameraGeometry* camera)
{
    double minX = 0;
    double maxX = 0;
    double minY = 0;
    double maxY = 0;
    std::unordered_map<double, int> xBuckets;
    std::unordered_map<double, int> yBuckets;

    for (auto sim : camera->recoData)
    {
        auto userInfo = static_cast<TTree*>(sim->Get("metadata"))->GetUserInfo();
        std::unordered_map<std::string, double> map;
        for (auto entry : *userInfo)
        {
            if (TString(entry->ClassName()) != "TParameter<double>") { continue; }
            auto e = static_cast<TParameter<double>*>(entry);
            map[std::string(e->GetName())] = e->GetVal();
        }
        if (map["sourcePosX"] > maxX) { maxX = map["sourcePosX"]; }
        if (map["sourcePosX"] < minX) { minX = map["sourcePosX"]; }
        if (map["sourcePosY"] > maxY) { maxY = map["sourcePosY"]; }
        if (map["sourcePosY"] < minY) { minY = map["sourcePosY"]; }

        // unsafe float comparison(may fail depending on simulation files)
        xBuckets[map["sourcePosX"]] += 1;
        yBuckets[map["sourcePosY"]] += 1;
    }

    double segSizeX = (maxX - minX) / (xBuckets.size() - 1);
    double segSizeY = (maxY - minY) / (yBuckets.size() - 1);

    camera->source.binX = xBuckets.size();
    camera->source.binY = yBuckets.size();
    camera->source.binZ = 1;
    camera->source.xRange = {minX - segSizeX / 2, maxX + segSizeX / 2};
    camera->source.yRange = {minY - segSizeY / 2, maxY + segSizeY / 2};
    camera->source.zRange = {-0.1, +0.1};
}
