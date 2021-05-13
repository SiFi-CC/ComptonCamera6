#include "CCMLEM.hh"
#include <CmdLineConfig.hh>
#include "IsectionPoint.hh"
#include <iostream>
using namespace std;

int main(int argc, char* argv[]) {
    
    CmdLineOption _path("output_path", "-op", "output path", "./results/");
    
    CmdLineConfig::instance()->ReadCmdLine(argc, argv);
    
    TString path = CmdLineOption::GetStringValue("output_path");
    
  /*  
    if (argc != 2) {
        cout << "To run type: ./mlem path_to_config" << endl;
        return 0;
    }*/
//   CmdLineOption opt_output("Output", "-o",
//                            "Output file (string), default: Hmatrix.root",
//                            "Hmatrix.root");

//  CmdLineConfig::instance()->ReadCmdLine(argc, argv);
//  TString outputfile(opt_output.GetStringValue()); '

    //TString path(argv[1]);
  
    CCMLEM* rec;

    try {
        
        rec = new CCMLEM(path);
        
    } catch (const char* message) {
        
        cout << message << endl;
        return 0;
    }

    rec->Reconstruct(1);
 // rec->HmatrixToFile(outputfile);

    delete rec;

    return 1;
}
