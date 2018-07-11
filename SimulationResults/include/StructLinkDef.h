#ifdef __CINT__

#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 

#pragma link C++ class std::vector<TVector3>+; 

#pragma link C++ struct PhysicVar+;
#pragma link C++ struct PhysicVec+;

#pragma link C++ class std::vector<PhysicVar>+; 
#pragma link C++ class std::vector<PhysicVec>+; 
#pragma link C++ class std::vector<PhysicVar*>+; 
#pragma link C++ class std::vector<PhysicVec*>+; 
#endif //__CINT__
