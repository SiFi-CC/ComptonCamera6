#include "InputReader.hh"
#include "InputReaderSimple.hh"

int main(){
 /*
  InputReader *in;
  
  try{
    in = new InputReader("../sources/results/results.root",false);
  }
  catch(const char *message){
   cout << message << endl;
   return 0;
  }
  
  in->AccessTree("tree_ft");
  in->Print(); 
  
  in->LoadEvent(1);*/
  //double x = in->GetValue();
  //for(int i=0; i<100; i++){
   //in->LoadEvent(i+1); 
  //}
  
  //TVector3 point1;
  //TVector3 point2;
  //double en1, en2;
  //in->ReadEvent(1,point1,point2,en1,en2);
  
  //point1.Print();
  //point2.Print();
  //cout << en1 << "\t" << en2 << endl;
  //-----

  InputReaderSimple *ins;
  
  try{  
    ins = new InputReaderSimple("../sources/results/CCSimulation_gen5.root");
  }
  catch(const char *message){
   cout << message << endl;
   return 0;
  }
  
  ins->Print();
  TVector3 *p0;
  
  for(int i=0; i<100; i++){
    ins->LoadEvent(i);
    p0 = ins->GetSourcePosition();
    p0->Print();
  }
  //-----
  //delete in;
  delete ins;
  
  return 1;
}
