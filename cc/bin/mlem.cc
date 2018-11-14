#include "CCMLEM.hh"
#include "IsectionPoint.hh"
#include <iostream>
using namespace std;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    cout << "To run type: ./mlem path_to_config" << endl;
    return 0;
  }

  TString path(argv[1]);

  try {
    CCMLEM rec(path);
    rec.Reconstruct();
  } catch (const char* message) {
    cout << message << endl;
    return 1;
  }

  return 0;
}
