#ifndef __DataStructConvert_H_
#define __DataStructConvert_H_ 1

#include "CLog.hh"
#include <TH2.h>
#include <TMatrixT.h>

namespace SiFi {
namespace tools {

TMatrixT<Double_t> convertHistogramToMatrix(TH2F* hist);
TH2F convertMatrixToHistogram(const char* name, const char* title,
                              TMatrixT<Double_t> matrix);
//
// Convert to vector by stacking columns one on the other
//
template <typename T> TMatrixT<T> vectorizeMatrix(const TMatrixT<T>& mat2D) {
  Int_t nRows = mat2D.GetNrows();
  Int_t nCols = mat2D.GetNcols();
  TMatrixT<Double_t> matVec(nRows * nCols, 1);

  for (int row = 0; row < nRows; row++) {
    for (int col = 0; col < nCols; col++) {
      matVec(row * nCols + col, 0) = mat2D(row, col);
    }
  }

  return matVec;
}

//
// Convert back to Matrix by assumong that columns are stacked one on another
template <typename T>
TMatrixT<T> unvectorizeMatrix(const TMatrixT<T>& matVec, Int_t nRows,
                              Int_t nCols) {
  if (matVec.GetNcols() != 1) {
    spdlog::error("unvectorizeMatrix: you need to pass column matrix");
    throw "wrong matrix dimensions";
  }
  if (matVec.GetNrows() != nRows * nCols) {
    spdlog::error("unvectorizeMatrix: number of columns needs to be equal to "
                  "nRows*nCols");
    throw "wrong matrix dimensions";
  }

  TMatrixT<T> mat2D(nRows, nCols);

  for (int row = 0; row < nRows; row++) {
    for (int col = 0; col < nCols; col++) {
      mat2D(row, col) = matVec(row * nCols + col, 0);
    }
  }

  return mat2D;
}

} // namespace tools
} // namespace SiFi

#endif
