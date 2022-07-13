#include "CLog.hh"

#include <TStopwatch.h>

#include "CmdLineConfig.hh"

#include "CMReconstruction.hh"
#include "G4SimulationAdapter.hh"

int main(int argc, char** argv)
{
    spdlog::set_level(spdlog::level::info);

    CmdLineOption opt_hmatrix("Hmatrix", "-hmat", "File with H matrix, default: Calculate", "");
    CmdLineOption opt_outputfile("output", "-o", "Output file", CmdLineArg::kString);

    CmdLineArg cmdarg_inputfile("input", "Input file", CmdLineArg::kString);
    CmdLineArg cmdarg_iter("iter", "N iterations", CmdLineArg::kInt);

    CmdLineConfig::instance()->ReadCmdLine(argc, argv);
    const Positional& args = CmdLineConfig::GetPositionalArguments();

    if (opt_hmatrix.GetStringValue())
    {
        spdlog::info("Hmatrix file: {}", opt_hmatrix.GetStringValue());
    }
    else
    {
        spdlog::info("Hmatrix will be calculated");
    }

    TString filename(args.at("input")->GetStringValue());
    Int_t iterations(args.at("iter")->GetIntValue());

    CMReconstruction reconstruction(filename);
    reconstruction.RunReconstruction(iterations);
    if (opt_outputfile.GetStringValue()) { reconstruction.Write(opt_outputfile.GetStringValue()); }
    else
    {
        reconstruction.Write(filename.ReplaceAll(".root", "_reconstruct.root"));
    }

    spdlog::info("Finished simulation");

    return 0;
}
