#include "CLog.hh"

#include <TFile.h>
#include <TString.h>

#include "CCMLEM.hh"
#include "CmdLineConfig.hh"


int main(int argc, char** argv) {
  //spdlog::set_level(spdlog::level::info);
  if (argc != 2) {
    cout << "To run type: ./smatrix path_to_config" << endl;
    return 0;
  }
  CmdLineOption opt_output("Output", "-o",
                           "Output file: Smatrix.root",
                           "Smatrix.root");

  CmdLineConfig::instance()->ReadCmdLine(argc, argv);

  
  TString outputfile(opt_output.GetStringValue());       // OUTPUT
  TString path(argv[1]);
  CCMLEM* reco;
  
  try {
    reco = new CCMLEM(path);
  } catch (const char* message) {
    cout << message << endl;
    return 0;
  }
  reco->Reconstruct(0);
  reco->SmatrixToFile(outputfile);
  
  delete reco;

  return 1;
} 
