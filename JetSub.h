#ifndef JETSUB
#define JETSUB

#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TreeStruct.h"
#include "TH2D.h"

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

class JetSub {
public:
  JetSub(const char* inputpath, const char* outputpath);
  ~JetSub();
  void Run();
  void Initialize();

private:
  TChain* inputfile;
  ExRootTreeReader* reader;
  TClonesArray* b_Jet4;
  TClonesArray* b_Jet6;
  TClonesArray* b_Jet10;
  TClonesArray* b_Particle;

  TFile* outputfile;
  TTree* writer;
  TClonesArray* o_Jet4;
  TClonesArray* o_Jet6;
  TClonesArray* o_Jet10;

  class QjetVar {
  public:
    QjetVar(): valid(false), mean(0), std(0) {}
    Bool_t valid;
    Float_t mean;
    Float_t std;
  };
  
  class NSubVar {
  public:
    NSubVar(Float_t t1 = 0, Float_t t2 = 0, Float_t t3 = 0): tau1(t1), tau2(t2), tau3(t3) {}
    Float_t tau1;
    Float_t tau2;
    Float_t tau3;
  };
  
  vector<int> pass_list4;
  vector<int> pass_list6;
  vector<int> pass_list10;
  
  vector<QjetVar> qvar4;
  vector<QjetVar> qvar6;
  vector<QjetVar> qvar10;

  vector<NSubVar> nvar4;
  vector<NSubVar> nvar6;
  vector<NSubVar> nvar10;

  Float_t zcut;
  Float_t dcut_fctr;
  Float_t exp_min;
  Float_t exp_max;
  Float_t rigidity;
  Float_t truncation_fctr;
  Int_t njets;
  Float_t r0;
  
  vector<int> Select(TClonesArray *jets);
  void Write();
  void WriteJet(TClonesArray* o_jet, TClonesArray* b_jet, vector<int>& pass_list, vector<QjetVar>& qvar, vector<NSubVar>& nvar);
  NSubVar Nsubjettiness(Jet* jet);
  QjetVar Qjet(Jet* jet);
};


#endif // JETSUB
