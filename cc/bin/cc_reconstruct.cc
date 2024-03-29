#include "CCReconstruction.hh"
#include <iostream>
using namespace std;

int main()
{

    Int_t nev = 100000;
    Int_t gen = 5;

    try
    {
        CCReconstruction rec(Form("../sources/results/CCSimulation_gen%i.root", gen),
                             Form("CCReconstruction_gen%i", gen), kFALSE);
        rec.RebuildSetupTxt();
        rec.ReconstructImage(0, nev);
    }
    catch (const char* message)
    {
        cout << message << endl;
        return 1;
    }

    return 0;
}
