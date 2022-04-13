#include "CCSimulation.hh"

#include "CLog.hh"

#include <CmdLineConfig.hh>

#include <TStopwatch.h>

#include <iostream>

#include <TMath.h>

using namespace std;

int main(int argc, char** argv)
{
    CmdLineOption _pos_y("pos_y", "-y", "y position of source", 0.0);
    CmdLineOption _pos_z("pos_z", "-z", "z position of source", 0.0);
    CmdLineOption _pos_x("pos_x", "-x", "x position of source", 0.0);
    CmdLineOption _number_events("no. of events", "-n", "number of events", 0);

    CmdLineOption _det_setup(
        "Setup", "-d", "Detector plane, 6 values required, default: 200:80:80:400:100:100", 0, 0);

    CmdLineOption _source_type("Source", "-s", "source type", -1000);
    CmdLineOption _angle_min("angle_min", "-amin", "minimum angle of generated photons", -TMath::Pi());
    CmdLineOption _angle_max("angle_max", "-amax", "maximum angle of generated photons", -TMath::Pi());
    CmdLineOption _output_path("OutputPath", "-opath", "output path", "./results/");

    CmdLineConfig::instance()->ReadCmdLine(argc, argv);

    Int_t nev = CmdLineOption::GetIntValue("no. of events");
    Double_t y = CmdLineOption::GetDoubleValue("pos_y");
    Double_t z = CmdLineOption::GetDoubleValue("pos_z");
    Double_t x = CmdLineOption::GetDoubleValue("pos_x");
    Double_t amin = CmdLineOption::GetDoubleValue("angle_min");
    Double_t amax = CmdLineOption::GetDoubleValue("angle_max");
    int gen = CmdLineOption::GetIntValue("Source");

    Double_t da = 200, db = 80, dc = 80, dd = 400, de = 100, df = 100;
    // Double_t da, db, dc, dd, de, df;
    if (_det_setup.GetArraySize() == 6)
    {
        da = _det_setup.GetDoubleArrayValue(1);
        db = _det_setup.GetDoubleArrayValue(2);
        dc = _det_setup.GetDoubleArrayValue(3);
        dd = _det_setup.GetDoubleArrayValue(4);
        de = _det_setup.GetDoubleArrayValue(5);
        df = _det_setup.GetDoubleArrayValue(6);
    }
    else if (_det_setup.GetArraySize() != 0)
    {
        spdlog::error("Detection setup - 6 parameters required, {} given",
                      _det_setup.GetArraySize());
        abort();
    }
    // DetPlane detector(da, db, dc, dd, de, df, "detector");
    TStopwatch t;
    t.Start();
    auto name = Form("CCSimulation_gen%i_corr_%.0f_%.0f_%.0f_no.%i_"
                     "scat._%.0fx%.0f_abs._%.0fx%.0f",
                     gen, x, y, z, nev, db, dc, de, df);
    printf("Creating CC simulation for '%s'\n", name);

    TString outputPath = CmdLineOption::GetStringValue("OutputPath");
    auto fGenVersion = CmdLineOption::GetIntValue("Source");

    auto sim = CCSimulation(name, outputPath, fGenVersion, amin, amax, kFALSE);

    // sim->SetGenVersion(gen.Atoi());
    sim.BuildSetup(da, db, dc, dd, de, df);
    sim.SetCoordinate(x, y, z);
    sim.Loop(nev);

    t.Stop();
    t.Print();

    return 0;
}
