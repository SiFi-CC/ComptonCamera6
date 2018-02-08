#include "PhysicsBase.hh"
#include "TMath.h"
#include "TRandom.h"
#include "TF1.h"
#include "TVector3.h"
#include <iostream>

using namespace std;

const double kR0  = 2.817E-15;	// m
const double kMe  = 0.510999;	// MeV/c2

//------------------------------------------------------------------
Double_t KleinNishina(Double_t *x, Double_t *par){
  
  //x[0]   - theta in deg
  //par[0] - energy of incident photon in MeV
  Double_t costheta = cos(x[0]*TMath::DegToRad());
  Double_t alpha    = par[0]/kMe;
  Double_t factor1  = (1. + costheta*costheta)/2.;
  Double_t factor2  = 1./(1. + alpha*(1. - costheta));
  Double_t factor3  = 1 + (alpha*alpha*(1 - costheta)*
                      (1 - costheta))/((1 + alpha*(1 - costheta))*
                      (1 + costheta*costheta));
 
  Double_t prob = kR0*kR0*factor1*factor2*factor2*factor3;	// m2
  //prob = prob*1E31;	// mb
 
  return prob;
}
//------------------------------------------------------------------

ClassImp(PhysicsBase);

//------------------------------------------------------------------
PhysicsBase::PhysicsBase(){
  SetName("physics");
  fTheta = 0.;
  fPhi = 0.;
}
//------------------------------------------------------------------
PhysicsBase::PhysicsBase(TString name){
  SetName(name);
  fTheta = 0.;
  fPhi = 0.;
}
//------------------------------------------------------------------
PhysicsBase::~PhysicsBase(){
  if(fFunction) delete fFunction;
}
//------------------------------------------------------------------
Double_t PhysicsBase::FindPhi(void){
  Double_t phi = gRandom->Uniform(-TMath::Pi(),TMath::Pi()); 	//rad 
  return phi;
}
//------------------------------------------------------------------
Double_t PhysicsBase::FindTheta(Double_t energy){
  fFunction = new TF1("fFunction",KleinNishina,0,180,1);
  fFunction->SetParameter(0,energy);
  Double_t tmp = fFunction->GetRandom(); 	//deg
  Double_t theta = tmp*TMath::DegToRad();	//rad
  return theta;
}
//------------------------------------------------------------------
Double_t PhysicsBase::NewEnergy(Double_t theta, Double_t initE){
  Double_t costheta = cos(theta);	//theta must be in rad
  Double_t alpha    = initE/kMe;
  Double_t finE     = initE/(1 + alpha*(1 - costheta));	//MeV
  return finE;
}
//------------------------------------------------------------------
Track* PhysicsBase::ComptonScatter(Track *initTrack, DetPlane *plane){
  
  TVector3 crossPoint;
  TVector3 finVersor;
  Bool_t   crossFlag;
  Double_t initE, finE;
  Double_t epsilon = 1.E-8;
  
  crossFlag = initTrack->FindCrossPoint(plane,crossPoint);
  if(crossFlag==kFALSE) return NULL;
  
  initE = initTrack->GetEnergy();
  fTheta = FindTheta(initE);	//rad
  //fTheta = TMath::Pi()/10.;
  fPhi = FindPhi();		//rad
  
  if(fabs(fTheta) < epsilon) 
    cout << "##### Warning! Theta angle after scattering still equals 0!" << endl;
  if(fabs(fPhi) < epsilon)
    cout << "##### Warning! Phi angle after scattering still equals 0!" << endl;
  
  finE = NewEnergy(fTheta,initE);
  
  //finVersor.SetXYZ(-1,0,0);
  //----- scattering
  TVector3 xPrim, yPrim, zPrim;
  TVector3 xVersor, yVersor;
  TVector3 xComp, yComp, zComp;
  xVersor.SetXYZ(1,0,0);
  yVersor.SetXYZ(0,1,0);
  zPrim = initTrack->GetVersor();
  yPrim = zPrim.Cross(yVersor);
  xPrim = yPrim.Cross(zPrim);
  yPrim.SetMag(1.);
  xPrim.SetMag(1.);
  xComp = cos(fPhi)*sqrt(1 - pow(cos(fTheta),2))*xPrim;
  yComp = sin(fPhi)*sqrt(1 - pow(cos(fTheta),2))*yPrim;
  zComp = cos(fTheta)*zPrim;
  finVersor = xComp + yComp + zComp;
  //-----
  
  Track *finTrack = new Track();
  finTrack->SetPoint(crossPoint);
  finTrack->SetEnergy(finE);
  finTrack->SetVersor(finVersor);
  finTrack->SetName("scattered");
  
  //----- theta check
  Double_t ang = initTrack->GetVersor().Angle(finVersor);
  if(fabs(ang-fTheta) > epsilon){
    cout << "##### Error in PhysicsBase::ComptonScatter!" << endl;
    cout << "##### Incorrect theta angle! PLease check!" << endl;
    cout << "Chosen theta = " << fTheta*TMath::RadToDeg() << " deg \t "
         << "Set theta = " << ang*TMath::RadToDeg() << " deg" << endl;
    return NULL;
  }
  //----- end of the theta check
  
  return finTrack;
}
//------------------------------------------------------------------
void PhysicsBase::Print(void){
 cout << "\nPhysicsBase::Print() for object " << GetName() << endl; 
 cout << "\tTheta scattering angle: " << fTheta*TMath::RadToDeg() << " deg" << endl;
 cout << "\tPhi scattering angle: " << fPhi*TMath::RadToDeg() << " deg" << endl; 
}
//------------------------------------------------------------------
