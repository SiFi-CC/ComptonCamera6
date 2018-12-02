#include "DR_GenerallStructs.hh"

PhysicVar::PhysicVar(double x, double y) : value(x), uncertainty(y) {}

PhysicVar PhysicVar::operator+(const PhysicVar& a) {
  value += a.value;
  uncertainty = sqrt(uncertainty * uncertainty + a.uncertainty * a.uncertainty);
  return *this;
}

PhysicVar PhysicVar::operator-(const PhysicVar& a) {
  value -= a.value;
  uncertainty = sqrt(uncertainty * uncertainty + a.uncertainty * a.uncertainty);
  return *this;
}

bool PhysicVar::operator==(const PhysicVar& a) const {
  return (value == a.value && uncertainty == a.uncertainty);
}

PhysicVar PhysicVar::Sum(const std::vector<PhysicVar> others) {
  for (int i = 0; i < others.size(); i++) {
    value += others[i].value;
    if (i == 0 && uncertainty != 0) {
      uncertainty = uncertainty * uncertainty +
                    others[i].uncertainty * others[i].uncertainty;
    } else {
      uncertainty += others[i].uncertainty * others[i].uncertainty;
    }
  }
  uncertainty = sqrt(uncertainty);
  return *this;
}

PhysicVar PhysicVar::Sum(const std::vector<PhysicVar*> others) {
  for (int i = 0; i < others.size(); i++) {
    value += others[i]->value;
    if (i == 0 && uncertainty != 0) {
      uncertainty = uncertainty * uncertainty +
                    others[i]->uncertainty * others[i]->uncertainty;
    } else {
      uncertainty += others[i]->uncertainty * others[i]->uncertainty;
    }
  }
  uncertainty = sqrt(uncertainty);
  return *this;
}

PhysicVec::PhysicVec(TVector3 x, TVector3 y) : position(x), uncertainty(y) {}

PhysicVec PhysicVec::operator+(const PhysicVec& a) {
  position += a.position;
  uncertainty.SetX(sqrt(uncertainty.x() * uncertainty.x() +
                        a.uncertainty.x() * a.uncertainty.x()));
  uncertainty.SetY(sqrt(uncertainty.y() * uncertainty.y() +
                        a.uncertainty.y() * a.uncertainty.y()));
  uncertainty.SetZ(sqrt(uncertainty.z() * uncertainty.z() +
                        a.uncertainty.z() * a.uncertainty.z()));
  return *this;
}
PhysicVec PhysicVec::operator-(const PhysicVec& a) {
  position -= a.position;
  uncertainty.SetX(sqrt(uncertainty.x() * uncertainty.x() +
                        a.uncertainty.x() * a.uncertainty.x()));
  uncertainty.SetY(sqrt(uncertainty.y() * uncertainty.y() +
                        a.uncertainty.y() * a.uncertainty.y()));
  uncertainty.SetZ(sqrt(uncertainty.z() * uncertainty.z() +
                        a.uncertainty.z() * a.uncertainty.z()));
  return *this;
}
bool PhysicVec::operator==(const PhysicVec& a) const {
  return (position == a.position && uncertainty == a.uncertainty);
}

PhysicVec PhysicVec::Sum(const std::vector<PhysicVec> others) {
  for (int i = 0; i < others.size(); i++) {
    position += others[i].position;
    if (i == 0 && uncertainty != TVector3(0, 0, 0)) {
      uncertainty.SetX(uncertainty.x() * uncertainty.x() +
                       others[i].uncertainty.x() * others[i].uncertainty.x());
      uncertainty.SetY(uncertainty.y() * uncertainty.y() +
                       others[i].uncertainty.y() * others[i].uncertainty.y());
      uncertainty.SetY(uncertainty.z() * uncertainty.z() +
                       others[i].uncertainty.z() * others[i].uncertainty.z());
    } else {
      uncertainty.SetX(uncertainty.x() +
                       others[i].uncertainty.x() * others[i].uncertainty.x());
      uncertainty.SetY(uncertainty.y() +
                       others[i].uncertainty.y() * others[i].uncertainty.y());
      uncertainty.SetY(uncertainty.z() +
                       others[i].uncertainty.z() * others[i].uncertainty.z());
    }
  }
  uncertainty.SetX(sqrt(uncertainty.x()));
  uncertainty.SetY(sqrt(uncertainty.y()));
  uncertainty.SetZ(sqrt(uncertainty.z()));
  return *this;
}

PhysicVec PhysicVec::Sum(const std::vector<PhysicVec*> others) {
  for (int i = 0; i < others.size(); i++) {
    position += others[i]->position;
    if (i == 0 && uncertainty != TVector3(0, 0, 0)) {
      uncertainty.SetX(uncertainty.x() * uncertainty.x() +
                       others[i]->uncertainty.x() * others[i]->uncertainty.x());
      uncertainty.SetY(uncertainty.y() * uncertainty.y() +
                       others[i]->uncertainty.y() * others[i]->uncertainty.y());
      uncertainty.SetY(uncertainty.z() * uncertainty.z() +
                       others[i]->uncertainty.z() * others[i]->uncertainty.z());
    } else {
      uncertainty.SetX(uncertainty.x() +
                       others[i]->uncertainty.x() * others[i]->uncertainty.x());
      uncertainty.SetY(uncertainty.y() +
                       others[i]->uncertainty.y() * others[i]->uncertainty.y());
      uncertainty.SetY(uncertainty.z() +
                       others[i]->uncertainty.z() * others[i]->uncertainty.z());
    }
  }
  uncertainty.SetX(sqrt(uncertainty.x()));
  uncertainty.SetY(sqrt(uncertainty.y()));
  uncertainty.SetZ(sqrt(uncertainty.z()));
  return *this;
}
