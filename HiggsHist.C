#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TImage.h"

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include "HiggsHist.h"
#include "HelperClasses.h"
#include "DelphesNTuple.h"

#include <string>
#include <iostream>

using namespace std;

class CHiggsHist : public CAnalysisData
{
private:
  TFile* file;
  TTree* minitTree;
  DelphesNTuple* mt;
  double xsection;

public:
  CHiggsHist();
  ~CHiggsHist();
  void Initialize(const char* inputFolder, const char* inputFile, const char* outputFolder, double xs);
  void ReInitialize(const char* inputFolder, const char* inputFile, double xs);
  void IterateOverEvents();
  void ProcessEvent();
};

CHiggsHist::CHiggsHist() {
	file = NULL;
	minitTree = NULL;
	mt = NULL;
}

CHiggsHist::~CHiggsHist() {
}

void CHiggsHist::Initialize(const char* inputFolder, const char* inputFile, const char* outputFolder, double xs) {
  xsection = xs;
  string inpFolder(inputFolder);
  string inpFile(inputFile);
  file = new TFile((inpFolder + inpFile).c_str());
  minitTree = (TTree*)file->Get("DelphesNTup");
  mt = new DelphesNTuple(minitTree);
  SimpleInitialize(inputFile, inputFolder, "res-");
}

void CHiggsHist::ReInitialize(const char* inputFolder, const char* inputFile, double xs) {
  xsection = xs;
  string inpFolder(inputFolder);
  file = new TFile((inpFolder + inputFile).c_str());
  minitTree = (TTree*)file->Get("DelphesNTup");
  mt = new DelphesNTuple(minitTree);
  // No need to call SimpleInitialize - the output file is already open
}

void CHiggsHist::IterateOverEvents() {
  for (Long64_t i = 0; i < mt->fChain->GetEntriesFast(); ++i) {
    Long64_t ientry = mt->LoadTree(i);
    if (ientry < 0) { 
      break;
    }
    
    mt->fChain->GetEntry(i);
    ProcessEvent();
  }
}

void CHiggsHist::ProcessEvent() {
  Fill("Flow Cut", 10);
  const double minBoostedJetPt = 2.0 * 126.0 / 1.0;
  // We store all events now, so check if there are exactly two leptons first
  if ((*mt->leptons_x).size() != 2) {
    return;
  }

  Fill("Flow Cut", 20);

  TLorentzVector l1((*mt->leptons_x)[0], (*mt->leptons_y)[0], (*mt->leptons_z)[0], (*mt->leptons_t)[0]);
  TLorentzVector l2((*mt->leptons_x)[1], (*mt->leptons_y)[1], (*mt->leptons_z)[1], (*mt->leptons_t)[1]);
  TLorentzVector dilepton = l1 + l2;

  Fill("m(ll)", dilepton.M());
  Fill("jets(r=0.4)", mt->jets_antikt_4_x->size());
  Fill("jets(r=0.6)", mt->jets_antikt_6_x->size());
  Fill("jets(r=1.0)", mt->jets_antikt_10_x->size());

  TLorentzVector boostedJet;
  TLorentzVector rjet1, rjet2;
  Fill("btag4", (double)(*mt->btags10)[0]);
  // If we are dealing with boosted jets
  if (mt->jets_antikt_4_x->size() < 2 || ((*mt->btags4)[0] & 2) == 0 || ((*mt->btags4)[1] & 2) == 0) {
    // have at least one boosted jet
    if (mt->jets_antikt_10_x->size() > 0 && ((*mt->btags10)[0] & 2) != 0) {
      boostedJet.SetXYZT((*mt->jets_antikt_10_x)[0], (*mt->jets_antikt_10_y)[0], (*mt->jets_antikt_10_z)[0], (*mt->jets_antikt_10_t)[0]);
      if (boostedJet.M() < 85.0 || boostedJet.M() > 165.0) {
	return;
      }

      Fill("Flow Cut", 40);
      if (boostedJet.Pt() < minBoostedJetPt) {
	return;
      }

      Fill("Flow Cut", 50);
      Fill("btag", (*mt->btags10)[0]);
      Fill("m(bb),b", boostedJet.M());
      Fill("pt(bb),b", boostedJet.Pt());
      Fill("m(llbb),b", (boostedJet + dilepton).M());
      Fill("deltaR(h,z),b", boostedJet.DeltaR(dilepton));
      Fill("deltaPhi(h,z),b", boostedJet.DeltaPhi(dilepton));
      Fill("pt(llbb),b", (boostedJet+dilepton).Pt());
    }
    
    Fill("Flow Cut", 30);
  } else {     
    rjet1.SetXYZT((*mt->jets_antikt_4_x)[0], (*mt->jets_antikt_4_y)[0], (*mt->jets_antikt_4_z)[0], (*mt->jets_antikt_4_t)[0]);
    rjet2.SetXYZT((*mt->jets_antikt_4_x)[1], (*mt->jets_antikt_4_y)[1], (*mt->jets_antikt_4_z)[1], (*mt->jets_antikt_4_t)[1]);
    TLorentzVector dijet = rjet2 + rjet1;
    if (dijet.M() < 85.0 || dijet.M() > 165.0) {
      return;
    }

    Fill("Flow Cut", 40);
    Fill("Flow Cut", 50);
    Fill("m(bb),r", dijet.M());
    Fill("m(llbb),r", (dijet + dilepton).M());
    Fill("deltaR(h,z),r", dijet.DeltaR(dilepton));
    Fill("deltaPhi(h,z),r", dijet.DeltaPhi(dilepton));
    Fill("pt(llbb),r", (dijet+dilepton).Pt());
    Fill("pt(bb),r", dijet.Pt());
  }
}

void HiggsHist(const char* inputFolder, std::vector<const char*> inputFiles, const char* outputFolder) {
  CHiggsHist hh;
  std::vector<const char*>::iterator i = inputFiles.begin();
  hh.Initialize(inputFolder, *i, outputFolder, 1.0);
  hh.IterateOverEvents();
  
  while (++i != inputFiles.end()) {
    hh.ReInitialize(inputFolder, *i, 1.0);
    hh.IterateOverEvents();
  }

  hh.SaveResults();
}

void HiggsHist(const char* inputFolder, const char* inputFile, const char* outputFolder) {
  CHiggsHist hh;
  hh.Initialize(inputFolder, inputFile, outputFolder, 1.0);
  hh.IterateOverEvents();
  hh.SaveResults();
}
