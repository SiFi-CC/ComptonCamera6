#include "MultiPointSource.hh"
#include "CLog.hh"
#include <TRandom.h>

// -----------------------
// MultiPointSource
// -----------------------

MultiPointSource::MultiPointSource(const TString fname)
{
    SetName("MultiPointSource");
    fInFileName = fname;
    if (!fInFileName.EndsWith(".mac"))
        throw "MultiPointSource constructor: not appropriate input file format, "
              ".mac file expected...";

    std::ifstream infile;
    infile.open(fInFileName.Data(), std::ios::in);
    if (!infile.is_open() || infile.fail() || infile.bad())
    {
        std::cout << "Input file " << fInFileName << " not opened correctly or corrupted"
                  << std::endl;
        std::cout << infile.is_open() << " " << infile.fail() << " " << infile.bad() << std::endl;
        throw "MultiPointSource constructor: cannot build source";
    }

    TString line = "";
    while (!infile.eof())
    {
        while (!(line.Contains("####") && line.Contains("Point")))
            line.ReadLine(infile);
        PointSource source(TVector3(0, 0, 0), 0, kFALSE);
        Bool_t initres = source.Init(infile);
        if (initres) fSources.push_back(source);
        fIntensities.push_back(source.Intensity());
    }
    fNumOfSources = fSources.size();
    if (NormalizeIntensities() == kFALSE)
        throw "MultiPointSource::NormalizeIntensities(): Sum of intensities of all "
              "sources = 0, check your input file...";
}

Track MultiPointSource::GenerateEvent()
{
    Int_t n = DiceSource();
    Track track = fSources[n].GenerateEvent();
    fhTheta->Fill(track.GetVersor().Theta());
    fhPhi->Fill(track.GetVersor().Phi());
    fhZY->Fill(track.GetPoint().Z(), track.GetPoint().Y());
    return track;
}

void MultiPointSource::AddSourceElement(const PointSource& source)
{
    fSources.push_back(source);
    fIterator = fSources.begin();
}

Bool_t MultiPointSource::NormalizeIntensities()
{
    Double_t sum = 0;
    for (Double_t f : fIntensities)
        sum = sum + f;
    if (sum < 1e-12) return kFALSE;
    for (int i = 0; i < fNumOfSources; i++)
        fIntensities[i] = fIntensities[i] / sum;
    return kTRUE;
}

Int_t MultiPointSource::DiceSource()
{
    Double_t rnm = gRandom->Uniform(0, 1);
    Double_t suml = 0;
    Double_t sumu = 0;
    for (int i = 0; i < fNumOfSources; i++)
    {
        sumu = sumu + fIntensities[i];
        if (rnm < sumu) return i;
    }
    return fNumOfSources;
}

void MultiPointSource::Print()
{
    std::cout << "This MultiPointSource contains " << fNumOfSources
              << " point-like sources:" << std::endl;
    for (int i = 0; i < fNumOfSources; i++)
    {
        std::cout << " \n\n### Source " << i + 1;
        fSources[i].Print();
    }
}
