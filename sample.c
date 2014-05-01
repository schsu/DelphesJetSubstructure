#include <TChain.h>
#include <TClonesArray.h>
#include <external/ExRootAnalysis/ExRootTreeReader.h>
#include "classes/DelphesClasses.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main() {
  TChain* chain = new TChain("Delphes");
  chain->Add("/Research/outputs/duplicates-10-fix.root");
  ExRootTreeReader* reader = new ExRootTreeReader(chain);
  TClonesArray* branchJet = reader->UseBranch("Jet10");
  TClonesArray* branchEFlowTrack = reader->UseBranch("EFlowTrack");
  TClonesArray* branchEFlowTower = reader->UseBranch("EFlowTower");
  TClonesArray* branchEFlowMuon = reader->UseBranch("EFlowMuon");
  int eventCount = reader->GetEntries();

  for (int i = 0; i < eventCount; ++i) {
    reader->ReadEntry(i);
    for (int jetNum = 0; jetNum < branchJet->GetEntries(); jetNum++) {
      Jet* jet = (Jet*)branchJet->At(jetNum);
      vector<TLorentzVector> constituents;
      for (int ci = 0; ci < jet->Constituents.GetEntriesFast(); ci++) {
	TObject* object = jet->Constituents.At(ci);
	if (object->IsA() == Muon::Class()) {
	  constituents.push_back(((Muon*)object)->P4());
	}
      }

      for (vector<TLorentzVector>::iterator it = constituents.begin(); it != constituents.end(); ++it) {
	for (vector<TLorentzVector>::iterator jt = it+1; jt != constituents.end(); ++jt) {
	  if(it->DeltaR(*jt) == 0.0) {
	    cout << "dupe" << endl;
	  }
	}
      }
    }
  }
  
  return 0;
}
