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
#include <map>
#include <algorithm>

typedef std::pair<TLorentzVector, int> jetInfo;

class CHiggsHist : public CAnalysisData
{
private:
  TFile* file;
  TTree* minitTree;
  DelphesNTuple* mt;
  double xsection;
  int eventCount;

  std::vector<TLorentzVector> ConvertBranchToVector(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& t);
  std::vector<jetInfo> ConvertBranchToJetVector(std::vector<double>& x, std::vector<double>& y, std::vector<double>&z, std::vector<double>& t, std::vector<int>& btag);
  std::vector<TLorentzVector> SelectGoodLeptons(std::vector<TLorentzVector>& leptons, double minPt, double maxEta);
  std::vector<jetInfo> SelectGoodJets(std::vector<jetInfo>& jets, std::vector<TLorentzVector>& electrons);
public:
  double fullEventCount;
  double amass;
  
  CHiggsHist();
  ~CHiggsHist();
  void Initialize(const char* inputFolder, const char* inputFile, const char* outputFolder, double xs);
  void ReInitialize(const char* inputFolder, const char* inputFile, double xs);
  void IterateOverEvents();
  void ProcessEvent();
  void ProcessHistograms();
  void DisplayFlow(const std::string& name, TH1F* hist);
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

std::vector<TLorentzVector>
CHiggsHist::SelectGoodLeptons(std::vector<TLorentzVector>& leptons, double minPt, double maxEta) {
  std::vector<TLorentzVector> goodLeptons;
  for (TLorentzVector& lv : leptons) {
    if (lv.Pt() > minPt && std::fabs(lv.Eta()) < maxEta) {
      goodLeptons.push_back(lv);
    }
  }

  return goodLeptons;
}

std::vector<jetInfo> 
CHiggsHist::SelectGoodJets(std::vector<jetInfo>& jets, std::vector<TLorentzVector>& electrons) {
  double minJetPt = 25.0;
  double maxJetEta = 2.0;

  std::vector<jetInfo> goodJets;
  for (jetInfo& j : jets) {
    if (j.first.Pt() < minJetPt || std::fabs(j.first.Eta()) > maxJetEta) {
      continue;
    }

    std::vector<TLorentzVector>::iterator match = std::find_if(electrons.begin(), electrons.end(), [j](const TLorentzVector& e) { return e.DeltaR(j.first) < 0.1; });
    if (match != electrons.end()) {
      // electrons are also detected as jets
      // if we have a match here, we don't want to treat it as a jet
      continue;
    }

    goodJets.push_back(j);
  }

  return goodJets;
}

std::vector<TLorentzVector>
CHiggsHist::ConvertBranchToVector(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector<double>& t) {
  std::vector<TLorentzVector> result;
  for (size_t i = 0; i < x.size(); ++i) {
    TLorentzVector lv(x[i], y[i], z[i], t[i]);
    result.push_back(lv);
  }

  return result;
}

std::vector<jetInfo> 
CHiggsHist::ConvertBranchToJetVector(std::vector<double>& x, std::vector<double>& y, std::vector<double>&z, std::vector<double>& t, std::vector<int>& btag) {
  std::vector<jetInfo> result;
  for (size_t i = 0; i < x.size(); ++i) {
    TLorentzVector lv(x[i], y[i], z[i], t[i]);
    result.push_back(jetInfo(lv, btag[i]));
  }

  return result;
}

void CHiggsHist::ProcessEvent() {
  double btagEfficiency = 0.6;
  double leptonMinPt = 5.0;
  double leptonMaxEta = 2.5;
  eventCount++;

  std::vector<TLorentzVector> electrons = ConvertBranchToVector(*mt->electrons_x, *mt->electrons_y, *mt->electrons_z, *mt->electrons_t);
  std::vector<TLorentzVector> muons = ConvertBranchToVector(*mt->muons_x, *mt->muons_y, *mt->muons_z, *mt->muons_t);

  std::vector<TLorentzVector> goodElectrons = SelectGoodLeptons(electrons, leptonMinPt, leptonMaxEta);
  std::vector<TLorentzVector> goodMuons = SelectGoodLeptons(muons, leptonMinPt, leptonMaxEta);

  Fill("Flow Cut", 10);

  bool matchingLeptons = (goodElectrons.size() == 2 && goodMuons.size() == 0) || (goodElectrons.size() == 0 && goodMuons.size() == 2);
  if (!matchingLeptons) {
    return;
  }

  Fill("Flow Cut", 20);

  std::vector<TLorentzVector> dilepton;
  if (goodElectrons.size() == 2) {
    dilepton = goodElectrons;
  } else {
    dilepton = goodMuons;
  }

  if (dilepton.size() != 2) {
    std::cout << "Something went wrong - there shoudl be exactly two leptons here." << std::endl;
    return;
  }

  double dileptonMass = (dilepton[0]+dilepton[1]).M();
  double dileptonPt = (dilepton[0]+dilepton[1]).Pt();
  if (dileptonMass < 80.0 || dileptonMass > 100.0) {
    return;
  }

  Fill("Flow Cut", 30);

  const double minBoostedJetPt = 2.0 * 126.0 / 1.2;

  bool resolved = false;
  std::vector<jetInfo> jetsLow = ConvertBranchToJetVector(*mt->jets_antikt_4_x, *mt->jets_antikt_4_y, *mt->jets_antikt_4_z, *mt->jets_antikt_4_t, *mt->btags4);
  std::vector<jetInfo> jetsMid = ConvertBranchToJetVector(*mt->jets_antikt_6_x, *mt->jets_antikt_6_y, *mt->jets_antikt_6_z, *mt->jets_antikt_6_t, *mt->btags6);
  std::vector<jetInfo> jetsHigh = ConvertBranchToJetVector(*mt->jets_antikt_10_x, *mt->jets_antikt_10_y, *mt->jets_antikt_10_z, *mt->jets_antikt_10_t, *mt->btags10);

  jetsLow = SelectGoodJets(jetsLow, goodElectrons);
  jetsMid = SelectGoodJets(jetsMid, goodElectrons);
  jetsHigh = SelectGoodJets(jetsHigh, goodElectrons);

  int btagFlag = 1;

  if (jetsLow.size() >= 2) {
    Fill("Flow Cut", 40);
    if (((jetsLow[0].second & btagFlag) != 0) && ((jetsLow[1].second & btagFlag) != 0)) {

      Fill("Flow Cut", 50);

      TLorentzVector dijet = jetsLow[0].first + jetsLow[1].first;
      if (dijet.M() > 126.0-15.0 && dijet.M() < 126.0+15.0) {
	Fill("Flow Cut", 60);
	Fill("m(ll),r", dileptonMass, xsection);
	Fill("pT(ll),r", dileptonPt, xsection);
	Fill("m(bb),r", dijet.M(), xsection);
	Fill("pT(bb),r", dijet.Pt(), xsection);
	Fill("m(llbb),r", (dijet+dilepton[0]+dilepton[1]).M(), xsection);
	Fill("pT(llbb),r", (dijet+dilepton[0]+dilepton[1]).Pt(), xsection);
	resolved = true;
      }
    }
  }

  if (!resolved) {
    Fill("Flow Cut", 70);
    std::vector<jetInfo> jetFilterPt;
    std::copy_if(jetsMid.begin(), jetsMid.end(), std::back_inserter(jetFilterPt), [](jetInfo& ji) { return ji.first.Pt() > 2.0 * 126.0 / .8; });
    if (jetFilterPt.size() > 0) {
      Fill("Flow Cut", 80);
      std::vector<jetInfo> jetFilterBTag;
      std::copy_if(jetFilterPt.begin(), jetFilterPt.end(), std::back_inserter(jetFilterBTag), [btagFlag](jetInfo& ji) { return ((ji.second & btagFlag) != 0); });
      if (jetFilterBTag.size() > 0) {
	Fill("Flow Cut", 90);
	if (jetFilterBTag[0].first.M() > 106.0 && jetFilterBTag[0].first.M() < 146.0) {
	  Fill("Flow Cut", 100);
	  Fill("m(ll),bl", dileptonMass, xsection);
	  Fill("pT(ll),bl", dileptonPt, xsection);
	  Fill("m(bb),bl", jetFilterBTag[0].first.M(), xsection);
	  Fill("pT(bb),bl", jetFilterBTag[0].first.Pt(), xsection);
	  Fill("m(llbb),bl", (jetFilterBTag[0].first+dilepton[0]+dilepton[1]).M(), xsection);
	  Fill("pT(llbb),bl", (jetFilterBTag[0].first+dilepton[0]+dilepton[1]).Pt(), xsection);
	}
      }
    }
  }

  if (!resolved) {
    Fill("Flow Cut", 110);
    std::vector<jetInfo> jetFilterPt;
    std::copy_if(jetsHigh.begin(), jetsHigh.end(), std::back_inserter(jetFilterPt), [](jetInfo& ji) { return ji.first.Pt() > 2.0 * 126.0 / 1.2; });
    if (jetFilterPt.size() > 0) {
      Fill("Flow Cut", 120);
      std::vector<jetInfo> jetFilterBTag;
      std::copy_if(jetFilterPt.begin(), jetFilterPt.end(), std::back_inserter(jetFilterBTag), [btagFlag](jetInfo& ji) { return ((ji.second & btagFlag) != 0); });
      if (jetFilterBTag.size() > 0) {
	Fill("Flow Cut", 130);
	if (jetFilterBTag[0].first.M() > 106.0 && jetFilterBTag[0].first.M() < 146.0) {
	  Fill("Flow Cut", 140);
	  Fill("m(ll),bh", dileptonMass, xsection);
	  Fill("pT(ll),bh", dileptonPt, xsection);
	  Fill("m(bb),bh", jetFilterBTag[0].first.M(), xsection);
	  Fill("pT(bb),bh", jetFilterBTag[0].first.Pt(), xsection);
	  Fill("m(llbb),bh", (jetFilterBTag[0].first+dilepton[0]+dilepton[1]).M(), xsection);
	  Fill("pT(llbb),bh", (jetFilterBTag[0].first+dilepton[0]+dilepton[1]).Pt(), xsection);
	}
      }
    }
  }

  return;
}

void
CHiggsHist::DisplayFlow(const std::string& name, TH1F* hist) {
  std::cout << "Flow " << name << std::endl;
  int nbins = hist->GetNbinsX();
  for (int i = 1; i <= nbins; ++i) {
    double items = hist->GetBinContent(i);
    if (items > 0.0) {
      std::cout << items << " ";
    }
  }

  std::cout << std::endl;
}

void
CHiggsHist::ProcessHistograms() {
  CAnalysisData::ProcessHistograms();
  for (auto& entry : histograms) {
    if (entry.first == "Flow Cut") {
      DisplayFlow("Resolved channel", entry.second);
    }
    
    entry.second->Scale(20000.0/(double)eventCount);
  #if 0
    if (entry.first == "pt(llbb),r") {
      gd.resolvedYield[amass] = entry.second->Integral();
    } else if (entry.first == "pt(llbb),bh") {
      gd.boostedHighYield[amass] = entry.second->Integral();
    } else if (entry.first == "pt(llbb),bl") {
      gd.boostedLowYield[amass] = entry.second->Integral();
    }
#endif
  }

  std::cout << eventCount << " events processed." << std::endl;
}

void HiggsHist(const char* inputFolder, std::vector<const char*>& inputFiles, std::vector<double>& crossSections, const char* outputFolder, double mass, double eventCountPerFile) {
  CHiggsHist hh;
  if (inputFiles.size() != crossSections.size()) {
    std::cout << "Input files number and cross sections number have to match!" << std::endl;
    return;
  }

  std::vector<const char*>::iterator i = inputFiles.begin();
  std::vector<double>::iterator ixs = crossSections.begin();
  hh.ReInitialize(inputFolder, *i, *ixs);
  hh.amass = mass;
  hh.fullEventCount = eventCountPerFile;
  hh.IterateOverEvents();
  
  while (++i != inputFiles.end()) {
    ++ixs;
    hh.fullEventCount += eventCountPerFile;
    hh.ReInitialize(inputFolder, *i, *ixs);
    hh.IterateOverEvents();
  }

  hh.SimpleInitialize("background.root", inputFolder, "res-");
  hh.SaveResults();
}

void HiggsHist(const char* inputFolder, const char* inputFile, const char* outputFolder, double xsection, double mass, double eventCount) {
  CHiggsHist hh;
  hh.Initialize(inputFolder, inputFile, outputFolder, xsection);
  hh.fullEventCount = eventCount;
  hh.amass = mass;
  hh.IterateOverEvents();
  hh.SaveResults();
}

