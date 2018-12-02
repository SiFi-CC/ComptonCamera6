#ifndef DR_STRUCTS_H
#define DR_STRUCTS_H

#include "TVector3.h"
#include <vector>
/*! Struct that builds a physical constant ( value + uncertainty)
 * If there is a PhysicVar example, value and uncertainty can be accessed via
 * example.value and example.uncertainty The operators for + and - are
 * overloaded so that gaussion error propagation is done
 * */
struct PhysicVar {
  double value;
  double uncertainty;

  PhysicVar(double x = 0, double y = 0);
  PhysicVar operator+(const PhysicVar& a);
  PhysicVar operator-(const PhysicVar& a);
  bool operator==(const PhysicVar& a) const;
  PhysicVar Sum(const std::vector<PhysicVar> others);
  PhysicVar Sum(const std::vector<PhysicVar*> others);
};
/*! Struct that builds a physical position (position + uncertainty)
 * If there is a PhysicVec example, value and uncertainty can be accessed via
 * example.value and example.uncertainty The operators for + and - are
 * overloaded so that gaussion error propagation is done
 */
struct PhysicVec {
  TVector3 position;
  TVector3 uncertainty;

  PhysicVec(TVector3 x = TVector3(0, 0, 0), TVector3 y = TVector3(0, 0, 0));
  PhysicVec operator+(const PhysicVec& a);
  PhysicVec operator-(const PhysicVec& a);
  bool operator==(const PhysicVec& a) const;
  PhysicVec Sum(const std::vector<PhysicVec> others);
  PhysicVec Sum(std::vector<PhysicVec*> others);
};

#endif // DR_STRUCTS_H
