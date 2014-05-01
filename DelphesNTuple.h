//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 29 13:58:33 2014 by ROOT version 5.34/18
// from TTree DelphesNTup/DelphesNTup
// found on file: /Research/outputs/mintree_jetsub_a-zh-triple-400GeV.root
//////////////////////////////////////////////////////////

#ifndef DelphesNTuple_h
#define DelphesNTuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.

class DelphesNTuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<double>  *leptons_x;
   vector<double>  *leptons_y;
   vector<double>  *leptons_z;
   vector<double>  *leptons_t;
   vector<double>  *jets_antikt_4_x;
   vector<double>  *jets_antikt_4_y;
   vector<double>  *jets_antikt_4_z;
   vector<double>  *jets_antikt_4_t;
   vector<double>  *jets_antikt_6_x;
   vector<double>  *jets_antikt_6_y;
   vector<double>  *jets_antikt_6_z;
   vector<double>  *jets_antikt_6_t;
   vector<double>  *jets_antikt_10_x;
   vector<double>  *jets_antikt_10_y;
   vector<double>  *jets_antikt_10_z;
   vector<double>  *jets_antikt_10_t;
   vector<double>  *tau1_4;
   vector<double>  *tau2_4;
   vector<double>  *tau3_4;
   vector<double>  *tau1_6;
   vector<double>  *tau2_6;
   vector<double>  *tau3_6;
   vector<double>  *tau1_10;
   vector<double>  *tau2_10;
   vector<double>  *tau3_10;
   vector<double>  *oldqj_m;
   vector<double>  *oldqj_m2;
   vector<double>  *oldqj_avg;
   vector<double>  *oldqj_rms;
   vector<double>  *oldqj_rmss;
   vector<double>  *newqj_m;
   vector<double>  *newqj_m2;
   vector<double>  *newqj_avg;
   vector<double>  *newqj_rms;
   vector<double>  *newqj_rmss;
   vector<int>     *btags4;
   vector<int>     *btags6;
   vector<int>     *btags10;

   // List of branches
   TBranch        *b_leptons_x;   //!
   TBranch        *b_leptons_y;   //!
   TBranch        *b_leptons_z;   //!
   TBranch        *b_leptons_t;   //!
   TBranch        *b_jets_antikt_4_x;   //!
   TBranch        *b_jets_antikt_4_y;   //!
   TBranch        *b_jets_antikt_4_z;   //!
   TBranch        *b_jets_antikt_4_t;   //!
   TBranch        *b_jets_antikt_6_x;   //!
   TBranch        *b_jets_antikt_6_y;   //!
   TBranch        *b_jets_antikt_6_z;   //!
   TBranch        *b_jets_antikt_6_t;   //!
   TBranch        *b_jets_antikt_10_x;   //!
   TBranch        *b_jets_antikt_10_y;   //!
   TBranch        *b_jets_antikt_10_z;   //!
   TBranch        *b_jets_antikt_10_t;   //!
   TBranch        *b_tau1_4;   //!
   TBranch        *b_tau2_4;   //!
   TBranch        *b_tau3_4;   //!
   TBranch        *b_tau1_6;   //!
   TBranch        *b_tau2_6;   //!
   TBranch        *b_tau3_6;   //!
   TBranch        *b_tau1_10;   //!
   TBranch        *b_tau2_10;   //!
   TBranch        *b_tau3_10;   //!
   TBranch        *b_oldqj_m;   //!
   TBranch        *b_oldqj_m2;   //!
   TBranch        *b_oldqj_avg;   //!
   TBranch        *b_oldqj_rms;   //!
   TBranch        *b_oldqj_rmss;   //!
   TBranch        *b_newqj_m;   //!
   TBranch        *b_newqj_m2;   //!
   TBranch        *b_newqj_avg;   //!
   TBranch        *b_newqj_rms;   //!
   TBranch        *b_newqj_rmss;   //!
   TBranch        *b_btags4;   //!
   TBranch        *b_btags6;   //!
   TBranch        *b_btags10;   //!

   DelphesNTuple(TTree *tree=0);
   virtual ~DelphesNTuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DelphesNTuple_cxx
DelphesNTuple::DelphesNTuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Research/outputs/mintree_jetsub_a-zh-triple-400GeV.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/Research/outputs/mintree_jetsub_a-zh-triple-400GeV.root");
      }
      f->GetObject("DelphesNTup",tree);

   }
   Init(tree);
}

DelphesNTuple::~DelphesNTuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DelphesNTuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DelphesNTuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DelphesNTuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   leptons_x = 0;
   leptons_y = 0;
   leptons_z = 0;
   leptons_t = 0;
   jets_antikt_4_x = 0;
   jets_antikt_4_y = 0;
   jets_antikt_4_z = 0;
   jets_antikt_4_t = 0;
   jets_antikt_6_x = 0;
   jets_antikt_6_y = 0;
   jets_antikt_6_z = 0;
   jets_antikt_6_t = 0;
   jets_antikt_10_x = 0;
   jets_antikt_10_y = 0;
   jets_antikt_10_z = 0;
   jets_antikt_10_t = 0;
   tau1_4 = 0;
   tau2_4 = 0;
   tau3_4 = 0;
   tau1_6 = 0;
   tau2_6 = 0;
   tau3_6 = 0;
   tau1_10 = 0;
   tau2_10 = 0;
   tau3_10 = 0;
   oldqj_m = 0;
   oldqj_m2 = 0;
   oldqj_avg = 0;
   oldqj_rms = 0;
   oldqj_rmss = 0;
   newqj_m = 0;
   newqj_m2 = 0;
   newqj_avg = 0;
   newqj_rms = 0;
   newqj_rmss = 0;
   btags4 = 0;
   btags6 = 0;
   btags10 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("leptons_x", &leptons_x, &b_leptons_x);
   fChain->SetBranchAddress("leptons_y", &leptons_y, &b_leptons_y);
   fChain->SetBranchAddress("leptons_z", &leptons_z, &b_leptons_z);
   fChain->SetBranchAddress("leptons_t", &leptons_t, &b_leptons_t);
   fChain->SetBranchAddress("jets_antikt_4_x", &jets_antikt_4_x, &b_jets_antikt_4_x);
   fChain->SetBranchAddress("jets_antikt_4_y", &jets_antikt_4_y, &b_jets_antikt_4_y);
   fChain->SetBranchAddress("jets_antikt_4_z", &jets_antikt_4_z, &b_jets_antikt_4_z);
   fChain->SetBranchAddress("jets_antikt_4_t", &jets_antikt_4_t, &b_jets_antikt_4_t);
   fChain->SetBranchAddress("jets_antikt_6_x", &jets_antikt_6_x, &b_jets_antikt_6_x);
   fChain->SetBranchAddress("jets_antikt_6_y", &jets_antikt_6_y, &b_jets_antikt_6_y);
   fChain->SetBranchAddress("jets_antikt_6_z", &jets_antikt_6_z, &b_jets_antikt_6_z);
   fChain->SetBranchAddress("jets_antikt_6_t", &jets_antikt_6_t, &b_jets_antikt_6_t);
   fChain->SetBranchAddress("jets_antikt_10_x", &jets_antikt_10_x, &b_jets_antikt_10_x);
   fChain->SetBranchAddress("jets_antikt_10_y", &jets_antikt_10_y, &b_jets_antikt_10_y);
   fChain->SetBranchAddress("jets_antikt_10_z", &jets_antikt_10_z, &b_jets_antikt_10_z);
   fChain->SetBranchAddress("jets_antikt_10_t", &jets_antikt_10_t, &b_jets_antikt_10_t);
   fChain->SetBranchAddress("tau1_4", &tau1_4, &b_tau1_4);
   fChain->SetBranchAddress("tau2_4", &tau2_4, &b_tau2_4);
   fChain->SetBranchAddress("tau3_4", &tau3_4, &b_tau3_4);
   fChain->SetBranchAddress("tau1_6", &tau1_6, &b_tau1_6);
   fChain->SetBranchAddress("tau2_6", &tau2_6, &b_tau2_6);
   fChain->SetBranchAddress("tau3_6", &tau3_6, &b_tau3_6);
   fChain->SetBranchAddress("tau1_10", &tau1_10, &b_tau1_10);
   fChain->SetBranchAddress("tau2_10", &tau2_10, &b_tau2_10);
   fChain->SetBranchAddress("tau3_10", &tau3_10, &b_tau3_10);
   fChain->SetBranchAddress("oldqj_m", &oldqj_m, &b_oldqj_m);
   fChain->SetBranchAddress("oldqj_m2", &oldqj_m2, &b_oldqj_m2);
   fChain->SetBranchAddress("oldqj_avg", &oldqj_avg, &b_oldqj_avg);
   fChain->SetBranchAddress("oldqj_rms", &oldqj_rms, &b_oldqj_rms);
   fChain->SetBranchAddress("oldqj_rmss", &oldqj_rmss, &b_oldqj_rmss);
   fChain->SetBranchAddress("newqj_m", &newqj_m, &b_newqj_m);
   fChain->SetBranchAddress("newqj_m2", &newqj_m2, &b_newqj_m2);
   fChain->SetBranchAddress("newqj_avg", &newqj_avg, &b_newqj_avg);
   fChain->SetBranchAddress("newqj_rms", &newqj_rms, &b_newqj_rms);
   fChain->SetBranchAddress("newqj_rmss", &newqj_rmss, &b_newqj_rmss);
   fChain->SetBranchAddress("btags4", &btags4, &b_btags4);
   fChain->SetBranchAddress("btags6", &btags6, &b_btags6);
   fChain->SetBranchAddress("btags10", &btags10, &b_btags10);
   Notify();
}

Bool_t DelphesNTuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DelphesNTuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DelphesNTuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DelphesNTuple_cxx
