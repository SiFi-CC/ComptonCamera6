#include "DataStructConvert.hh"

namespace SiFi {
namespace tools {

TMatrixT<Double_t> convertHistogramToMatrix(TH2F* hist) {
  Int_t nBinsX = hist->GetXaxis()->GetNbins();
  Int_t nBinsY = hist->GetYaxis()->GetNbins();
  TMatrixT<Double_t> matrix(nBinsY, nBinsX);

  for (int binX = 1; binX <= nBinsX; binX++) {
    for (int binY = 1; binY <= nBinsY; binY++) {
      // matrix(nBinsY - binY, binX - 1) = hist->GetBinContent(binX, binY);
      matrix(binY - 1, binX - 1) = hist->GetBinContent(binX, binY);
    }
  }
  return matrix;
}

TH2F convertMatrixToHistogram(const char* name, const char* title,
                              TMatrixT<Double_t> matrix,
                              std::pair<double, double> xRange,
                              std::pair<double, double> yRange) {
  Int_t nBinsX = matrix.GetNcols();
  Int_t nBinsY = matrix.GetNrows();
  TH2F hist(name, title, nBinsX, xRange.first, xRange.second, nBinsY,
            yRange.first, yRange.second);

  for (int binX = 1; binX <= nBinsX; binX++) {
    for (int binY = 1; binY <= nBinsY; binY++) {
      // hist.SetBinContent(binX, binY, matrix(nBinsY - binY, binX - 1));
      hist.SetBinContent(binX, binY, matrix(binY - 1, binX - 1));
    }
  }

  return hist;
}

} // namespace tools
} // namespace SiFi
