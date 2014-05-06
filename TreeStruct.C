#include "TLorentzVector.h"
#include "TreeStruct.h"

using namespace treestruct;

TLorentzVector Jet::P4() {
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}
