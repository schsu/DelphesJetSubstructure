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
  int eventCount;

public:
  CHiggsHist();
  ~CHiggsHist();
  void Initialize(const char* inputFolder, const char* inputFile, const char* outputFolder, double xs);
  void ReInitialize(const char* inputFolder, const char* inputFile, double xs);
  void IterateOverEvents();
  void ProcessEvent();
  void ProcessHistograms();
};

CHiggsHist::CHiggsHist() {
	file = NULL;
	minitTree = NULL;
	mt = NULL;
	eventCount = 0;
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
  std::cout << inputFile << std::endl;
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
  double btagEfficiency = 0.6;
  eventCount++;

  Fill("Flow Cut, r", 10);
  Fill("Flow Cut, b", 10);
  const double minBoostedJetPt = 2.0 * 126.0 / 1.0;
  // We store all events now, so check if there are exactly two leptons first
  if ((*mt->leptons_x).size() != 2) {
    return;
  }

  Fill("Flow Cut, r", 20);
  Fill("Flow Cut, b", 20);

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
    Fill("Flow Cut, b", 30);
    if (mt->jets_antikt_10_x->size() > 0 && ((*mt->btags10)[0] & 2) != 0) {
      boostedJet.SetXYZT((*mt->jets_antikt_10_x)[0], (*mt->jets_antikt_10_y)[0], (*mt->jets_antikt_10_z)[0], (*mt->jets_antikt_10_t)[0]);
      if (boostedJet.M() < 85.0 || boostedJet.M() > 165.0) {
	return;
      }

      Fill("Flow Cut, b", 40);
      if (boostedJet.Pt() < minBoostedJetPt) {
	return;
      }

      Fill("Flow Cut, b", 50);
      Fill("m(bb),b", boostedJet.M(), btagEfficiency * xsection);
      Fill("pt(bb),b", boostedJet.Pt(), btagEfficiency * xsection);
      Fill("m(llbb),b", (boostedJet + dilepton).M(), btagEfficiency * xsection);
      Fill("deltaR(h,z),b", boostedJet.DeltaR(dilepton), btagEfficiency * xsection);
      Fill("deltaPhi(h,z),b", boostedJet.DeltaPhi(dilepton), btagEfficiency * xsection);
      Fill("pt(llbb),b", (boostedJet+dilepton).Pt(), btagEfficiency * xsection);
    }
  } else {
    Fill("Flow Cut, r", 30);
    rjet1.SetXYZT((*mt->jets_antikt_4_x)[0], (*mt->jets_antikt_4_y)[0], (*mt->jets_antikt_4_z)[0], (*mt->jets_antikt_4_t)[0]);
    rjet2.SetXYZT((*mt->jets_antikt_4_x)[1], (*mt->jets_antikt_4_y)[1], (*mt->jets_antikt_4_z)[1], (*mt->jets_antikt_4_t)[1]);
    TLorentzVector dijet = rjet2 + rjet1;
    if (dijet.M() < 85.0 || dijet.M() > 165.0) {
      return;
    }

    Fill("Flow Cut, r", 40);
    Fill("m(bb),r", dijet.M(), btagEfficiency * btagEfficiency * xsection);
    Fill("m(llbb),r", (dijet + dilepton).M(), btagEfficiency * btagEfficiency * xsection);
    Fill("deltaR(h,z),r", dijet.DeltaR(dilepton), btagEfficiency * btagEfficiency * xsection);
    Fill("deltaPhi(h,z),r", dijet.DeltaPhi(dilepton), btagEfficiency * btagEfficiency * xsection);
    Fill("pt(llbb),r", (dijet+dilepton).Pt(), btagEfficiency * btagEfficiency * xsection);
    Fill("pt(bb),r", dijet.Pt(), btagEfficiency * btagEfficiency * xsection);
  }
}

void
CHiggsHist::ProcessHistograms() {
  CAnalysisData::ProcessHistograms();
  for (auto& entry : histograms) {
    entry.second->Scale(20000.0/(double)eventCount);
  }

  std::cout << eventCount << " events processed." << std::endl;
}

void HiggsHist(const char* inputFolder, std::vector<const char*>& inputFiles, std::vector<double>& crossSections, const char* outputFolder) {
  CHiggsHist hh;
  if (inputFiles.size() != crossSections.size()) {
    std::cout << "Input files number and cross sections number have to match!" << std::endl;
    return;
  }

  std::vector<const char*>::iterator i = inputFiles.begin();
  std::vector<double>::iterator ixs = crossSections.begin();
  hh.ReInitialize(inputFolder, *i, *ixs);
  hh.IterateOverEvents();
  
  while (++i != inputFiles.end()) {
    ++ixs;
    hh.ReInitialize(inputFolder, *i, *ixs);
    hh.IterateOverEvents();
  }

  hh.SimpleInitialize("background.root", inputFolder, "res-");
  hh.SaveResults();
}

void HiggsHist(const char* inputFolder, const char* inputFile, const char* outputFolder) {
  CHiggsHist hh;
  hh.Initialize(inputFolder, inputFile, outputFolder, 1.0);
  hh.IterateOverEvents();
  hh.SaveResults();
}

