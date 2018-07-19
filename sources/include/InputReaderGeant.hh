#ifndef __InputReaderGeant_H_
#define __InputReaderGeant_H_ 1
#include "InputReader.hh"
#include <iostream>

using namespace std;

class InputReaderGeant : public InputReader{
  
public:
  InputReaderGeant();
  InputReaderGeant(TString path);
  ~InputReaderGeant();
  
private:
  bool AccessTree(TString name);
  
  ClassDef(InputReaderGeant,0)
};

#endif
