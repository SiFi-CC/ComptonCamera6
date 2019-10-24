#include "CLog.hh"
#include <TFile.h>
#include <TString.h>
#include "CMReconstruction.hh"
#include "CmdLineConfig.hh"



int main(int argc, char** argv) {
	spdlog::set_level(spdlog::level::info);

	CmdLineOption opt_output("Output", "-o",
                           "Output file (string), default: Hmatrix.root",
                           "Hmatrix.root");

	CmdLineConfig::instance()->ReadCmdLine(argc, argv);
	PositionalArgs args = CmdLineOption::GetPositionalArguments();

	spdlog::info("Hello, Vitalii!");
	if (args.size() != 1 ) {
		spdlog::info(
		    "type: './cm_reconstruct [INPUTFILENAME] ' to start:\n\n"
		    "where:\n\n"
		    "INPUTFILENAME - is an input file from simulations\n\n");
		return 1;
	}		

	spdlog::info("Outputfile: {}",opt_output.GetStringValue());

    TString inputfile(args[0]); //INPUT
    TString outputfile(opt_output.GetStringValue()); //OUTPUT

	CMReconstruction reconstruction(inputfile);
	
	reconstruction.HmatrixToFile(outputfile,0);

}