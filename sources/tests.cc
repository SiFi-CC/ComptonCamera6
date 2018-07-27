#include "InputReader.hh"
#include "InputReaderSimple.hh"
#include "InputReaderGeant.hh"

int main(){
  
  TString fname = "../sources/results/CCSimulation_gen5.root";
  int istart = 10;
  int istop  = 100;
  
  cout << "===== Test of InputReader class" << endl;
  
  InputReader *in;
  
  try{
    in = new InputReader(fname);
  }
  catch(const char *message){
   cout << message << endl;
   return 0;
  }

  in->Print();
  
  //-----

  cout << "===== Test of InputReaderSimple class" << endl;
  
  InputReaderSimple *ins;
  
  try{  
    ins = new InputReaderSimple(fname);
  }
  catch(const char *message){
   cout << message << endl;
   return 0;
  }
  
  ins->Print();
  
  TVector3 *point0;
  TVector3 *point1;
  TVector3 *point2;
  TVector3 *versor1;
  TVector3 *versor2;
  double   energy0;
  double   energy1;
  double   energy2;
  
  for(int i=istart; i<istop; i++){
    ins->LoadEvent(i);
    point0 = ins->GetPositionPrimary();
    point1 = ins->GetPositionScattering();
    point2 = ins->GetPositionAbsorption();
    versor1 = ins->GetGammaDirPrimary();
    versor2 = ins->GetGammaDirScattered();
    energy0 = ins->GetEnergyPrimary();
    energy1 = ins->GetEnergyLoss();
    energy2 = ins->GetEnergyScattered();
  }
  
  //-----
  
  cout << "===== Test of InputReaderGeant class" << endl;
  
  fname = "../SimulationResults/test.root";
  InputReaderGeant *ing;
  
  try{
    ing = new InputReaderGeant(fname); 
  }
  catch(const char *message){
    cout << message << endl;
    return 0;
  }
  
  ing->Print();
  
  //-----
  delete in;
  delete ins;

  return 1;
}
