#ifndef TREESTRUCT
#define TREESTRUCT

#include "TObject.h"
#include "TLorentzVector.h"

namespace treestruct {
  class Jet: public TObject {
  public:
    Float_t PT;
    Float_t Eta;
    Float_t Phi;
    Float_t Mass;
    
    Float_t Tau1;
    Float_t Tau2;
    Float_t Tau3;
    
    Float_t Mean;
    Float_t Std;
    Float_t Volality;
    
    TLorentzVector P4();
    
    ClassDef(Jet, 1);
  };
  
};

#endif // TREESTRUCT
