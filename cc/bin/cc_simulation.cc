#include "CCSimulation.hh"
#include <TStopwatch.h>
#include <iostream>
#include <CmdLineConfig.hh>

using namespace std;

int main(int argc, char** argv){
  
  CmdLineOption _pos_y("pos_y", "-y", "y position of source" , 0.0);
  CmdLineOption _pos_z("pos_z", "-z", "z position of source" , 0.0);
  CmdLineOption _pos_x("pos_x", "-x", "x position of source" , 0.0);
  CmdLineOption _number_events("no. of events", "-n", "number of events" , 0);
  CmdLineConfig::instance()->ReadCmdLine(argc, argv);
  
//   if (argc != 2){
//     cout<< " To run type: './bin/cc_simulation fGenVersion' "<<endl;
//     return 0;
//   }
  
//   TString gen(argv[1]);

  //Int_t nev = 100;
  Int_t nev = CmdLineOption::GetIntValue("no. of events");
  Double_t y = CmdLineOption::GetDoubleValue("pos_y");
  Double_t z = CmdLineOption::GetDoubleValue("pos_z");
  Double_t x = CmdLineOption::GetDoubleValue("pos_x");
  int gen = CmdLineOption::GetIntValue("Source");

  CCSimulation *sim;  
  TStopwatch t;
  t.Start();  
//   for(int i=0; i<n; i++){
      
//       for(int j=0; j<n; j++){
	sim = new CCSimulation(
	  Form("CCSimulation_gen%i_corr_%.0f_%.0f_%.0f_no.%i",gen,x,y,z,nev),kFALSE);
// 	sim->SetGenVersion(gen.Atoi());
	sim->BuildSetup(200,80,80,400,100,100);
        sim->SetCoordinate(x,y,z);
        sim->Loop(nev);
	delete sim;
//       }
      
//     }
    
    t.Stop(); 
    t.Print();
 
  return 1; 
}

