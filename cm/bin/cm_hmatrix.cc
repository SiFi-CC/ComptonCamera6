#include "CLog.hh"

#include <TFile.h>
#include <TString.h>

#include "CMReconstruction.hh"
#include "CmdLineConfig.hh"

int main(int argc, char** argv)
{
    spdlog::set_level(spdlog::level::info);

    CmdLineOption opt_output("Output", "-o", "Output file (string), default: Hmatrix.root",
                             "Hmatrix.root");

    CmdLineOption opt_events("Events", "-n",
                             "Number of events for each vertex, default: 100000 (integer)", 100000);

    CmdLineArg cmdarg_input("input", "Input file", CmdLineArg::kString);

    CmdLineConfig::instance()->ReadCmdLine(argc, argv);

    const Positional& args = CmdLineConfig::GetPositionalArguments();

    spdlog::info("Outputfile: {}", opt_output.GetStringValue());

    TString inputfile(args.at("input")->GetStringValue()); // INPUT
    TString outputfile(opt_output.GetStringValue());       // OUTPUT

    CMReconstruction reconstruction(inputfile);

    reconstruction.HmatrixToFile(outputfile);
}