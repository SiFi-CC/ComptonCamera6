#include "SMElement.hh"
#include <iostream>

using namespace std;

ClassImp(SMElement);

//------------------------------------------------------------------
/// Default constructor.
SMElement::SMElement()
{
    feventno = 0;
    fGlobalBin = 0;
    fdist = 0;
}
//------------------------------------------------------------------
/// Standard constructor
///\param eventno (Int_t) - event number.
///\param bin (Int_t) - global bin number.
///\param dist (Double_t) - lenght of track.
SMElement::SMElement(Int_t eventno, Int_t bin, Double_t dist)
{
    feventno = eventno;
    fGlobalBin = bin;
    fdist = dist;
}
//------------------------------------------------------------------
/// Deafault destructor
SMElement::~SMElement() {}
//------------------------------------------------------------------
/// Sets event number and global bin number and length of track.
///\param eventno (Int_t) - event number.
///\param bin (Int_t) - global bin number.
///\param dist (Double_t) - length of track.
void SMElement::SetEvBinDist(Int_t eventno, Int_t bin, Double_t dist)
{
    feventno = eventno;
    fGlobalBin = bin;
    fdist = dist;
}
//------------------------------------------------------------------
/// Sets event number.
///\param i (Int_t) - event number.
void SMElement::SetEvent(Int_t i) { feventno = i; }
//------------------------------------------------------------------
/// Sets global bin number.
///\param b (Int_t) - global bin number.
void SMElement::SetBin(Int_t b) { fGlobalBin = b; }
//------------------------------------------------------------------
/// Sets length of track.
///\param d (Double_t) - length of track.
void SMElement::SetDist(Double_t d) { fdist = d; }
//------------------------------------------------------------------
/// Returns event number.
Int_t SMElement::GetEvent(void) const { return feventno; }
//------------------------------------------------------------------
/// Returns global bin number.
Int_t SMElement::GetBin(void) const { return fGlobalBin; }
//------------------------------------------------------------------
/// Returns lenght of track.
Double_t SMElement::GetDist(void) { return fdist; }
//------------------------------------------------------------------
/// Prints details of the SMElement class object.
void SMElement::Print(void)
{
    cout << "\nSMElement::Print()" << endl;
    cout << "\tEvent number: \t" << feventno << endl;
    cout << "\tGlobal bin number: \t" << fGlobalBin << endl;
    cout << "\tLength of track: " << fdist << endl;
}
//------------------------------------------------------------------
