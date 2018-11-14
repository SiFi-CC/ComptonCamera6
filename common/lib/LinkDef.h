 #ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Track;
#pragma link C++ class DetPlane;
#pragma link C++ class IsectionPoint;
#pragma link C++ class InputReader;
#pragma link C++ class InputReaderSimple;
#pragma link C++ class InputReaderGeant;

#pragma link C++ struct PhysicVar+;
#pragma link C++ struct PhysicVec+;

#pragma link C++ class std::vector<PhysicVar>+;
#pragma link C++ class std::vector<PhysicVec>+;
#pragma link C++ class std::vector<PhysicVar*>+;
#pragma link C++ class std::vector<PhysicVec*>+;

#endif