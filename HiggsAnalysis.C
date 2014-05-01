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

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <numeric>

#include "HelperClasses.h"

class CHiggsAnalysis : public CAnalysisData {
public:
  float scaleFactor;
  CHiggsAnalysis();
  virtual ~CHiggsAnalysis();

  virtual void ProcessEvent();
  int MaxParticleForJet(BaseParticle* jet, bool& hasBq, bool& hasHiggs);
  bool HasHiggsParent(BaseParticle* jet);
  virtual void ProcessHistograms();
  int trueBJet;
  int taggedBJet;
  void DoCounting();
  void SetSpecialHistogramRanges();

  void SelectGoodJets(std::vector<BaseParticle*>& goodJets, TClonesArray* branchJets, const std::vector<BaseParticle*>& goodElectrons);
  void ResetGoodArrays();

  float jetMinPt;
  float jetMaxEta;
  float minMassCut;
  float maxMassCut;
  bool useMassCut;
};

class CHiggsTruth : public CAnalysisData {
public:
  CHiggsTruth();
  virtual ~CHiggsTruth();
  float heavyHiggsMass;

  virtual void ProcessEvent();
};

CHiggsTruth::CHiggsTruth() {
}

CHiggsTruth::~CHiggsTruth() {
}

void CHiggsTruth::ProcessEvent() {
  int ahiggs = 0;

  for (int i = 2; i < branchParticle->GetEntries(); i++) {
    GenParticle* particle = (GenParticle*)branchParticle->At(i);
    if (particle->Status == 3 && particle->PID == 36) {
      ahiggs = i;
      break;
    }
  }

  int hhiggs = 0;
  int z = 0;
  vector<int> leptons;
  vector<int> bqs;
  if (ahiggs > 0) {
    // Find the sm higgs
    for (int i = 0; i < branchParticle->GetEntries(); i++) {
      GenParticle* particle = (GenParticle*)branchParticle->At(i);
      if (particle->PID == 25 && particle->M1 == ahiggs) {
        hhiggs = i;
        break;
      }
    }
    // Find the z
    for (int i = 0; i < branchParticle->GetEntries(); i++) {
      GenParticle* particle = (GenParticle*)branchParticle->At(i);
      if (particle->PID == 23 && particle->M1 == ahiggs) {
        z = i;
        break;
      }
    }

    if (hhiggs > 0 && z > 0) {
      for (int i = 0; i < branchParticle->GetEntries(); i++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(i);
        if ((abs(particle->PID) == 13 || abs(particle->PID) == 11) && particle->M1 == z) {
          leptons.push_back(i);
          if (leptons.size() == 2) {
            break;
          }
        }
      }

      for (int i = 0; i < branchParticle->GetEntries(); i++) {
        GenParticle* particle = (GenParticle*)branchParticle->At(i);
        if ((abs(particle->PID) == 5) && particle->M1 == hhiggs) {
          bqs.push_back(i);
          if (bqs.size() == 2) {
            break;
          }
        }
      }
    }

    if (leptons.size() < 2 || bqs.size() < 2) {
      cout << "bad" << endl;
    }

    TLorentzVector dilepton;
    TLorentzVector dijet;
    dilepton = GetParticle(leptons[0])->P4() + GetParticle(leptons[1])->P4();
    dijet = GetParticle(bqs[0])->P4() + GetParticle(bqs[1])->P4();
    TLorentzVector h;
    h = GetParticle(hhiggs)->P4();
    TLorentzVector a;
    a = GetParticle(ahiggs)->P4();
    TLorentzVector zv;
    zv = GetParticle(z)->P4();

    SetHistogramMinMax("m(bb)", 124.5, 125.5);
    SetHistogramMinMax("m(ll)", 50.0, 130.0);
    SetHistogramMinMax("m(llbb)", heavyHiggsMass - 2.5, heavyHiggsMass + 2.5);
    SetHistogramMinMax("m(h)", 124.5, 125.5);
    SetHistogramMinMax("m(a)", heavyHiggsMass - 2.5, heavyHiggsMass + 2.5);
    SetHistogramMinMax("m(z)", 50.0, 130.0);

    SetHistogramMinMax("pt(llbb)", 0.0, 1000.0);
    SetHistogramMinMax("pt(ll)", 0.0, 1000.0);
    SetHistogramMinMax("pt(bb)", 0.0, 1000.0);
    SetHistogramMinMax("pt(h)", 0.0, 1000.0);
    SetHistogramMinMax("pt(z)", 0.0, 1000.0);
    SetHistogramMinMax("pt(a)", 0.0, 1000.0);
    SetHistogramMinMax("pt(h)", 0.0, 1000.0);

    SetHistogramMinMax("eta(ll)", -6.0, 6.0);
    SetHistogramMinMax("eta(bb)", -6.0, 6.0);
    SetHistogramMinMax("eta(llbb)", -6.0, 6.0);
    SetHistogramMinMax("eta(h)", -6.0, 6.0);
    SetHistogramMinMax("eta(z)", -6.0, 6.0);
    SetHistogramMinMax("eta(a)", -6.0, 6.0);

    float normalization = 1.0 / 50000.0;

    Fill("m(ll)", dilepton.M(), normalization);
    Fill("pt(ll)", dilepton.Pt(), normalization);
    Fill("eta(ll)", dilepton.Eta(), normalization);

    Fill("m(a)", a.M(), normalization);
    Fill("pt(a)", a.Pt(), normalization);
    Fill("eta(a)", a.Eta(), normalization);

    Fill("m(h)", h.M(), normalization);
    Fill("pt(h)", h.Pt(), normalization);
    Fill("eta(h)", h.Eta(), normalization);

    Fill("m(z)", zv.M(), normalization);
    Fill("pt(z)", zv.Pt(), normalization);
    Fill("eta(z)", zv.Eta(), normalization);

    Fill("m(bb)", dijet.M(), normalization);
    Fill("pt(bb)", dijet.Pt(), normalization);
    Fill("eta(bb)", dijet.Eta(), normalization);

    TLorentzVector ah = dijet + dilepton;
    Fill("m(llbb)", ah.M(), normalization);
    Fill("pt(llbb)", ah.Pt(), normalization);
    Fill("eta(llbb)", ah.Eta(), normalization);
  }
}

class CBackGroundAnalysis : public CHiggsAnalysis {
private:
  vector<TString> inputFileList;
  vector<float> xsectionsList;
  TString currentInputFolder;
  size_t current;
  bool loadNextFile();
public:
  CBackGroundAnalysis();
  virtual ~CBackGroundAnalysis();
  void Initialize(const char* inputFolder, vector<TString>& inputFile, vector<float>& xsections, const char* outputFolder, const char* outputPrefix);
  virtual void IterateOverEvents();
};

class CHiggsAnalysisNoCuts : public CAnalysisData {
public:
  CHiggsAnalysisNoCuts();
  virtual ~CHiggsAnalysisNoCuts();

  virtual void ProcessEvent();
};

CHiggsAnalysis::CHiggsAnalysis() {
  trueBJet = 0;
  taggedBJet = 0;
  jetMinPt = 20.0;
  jetMaxEta = 2.5;
  scaleFactor = 50000.0;
}

CHiggsAnalysis::~CHiggsAnalysis() {
  cout << "True:" << trueBJet << " Tagged: " << taggedBJet << endl;
}

void CHiggsAnalysis::SetSpecialHistogramRanges() {
  SetHistogramMinMax("m(llbb) boosted", 0.0, 2000.0);
  SetHistogramMinMax("m(llbb) resolved", 0.0, 2000.0);
}

void CHiggsAnalysis::DoCounting() {
  for (int i = 0; i < branchJet4->GetEntries(); i++) {
    Jet* jet = (Jet*)GetJet4(i)->underlying();
    if (jet->BTag & 1) {
      taggedBJet++;
    }

    if (jet->BTag & 2) {
      trueBJet++;
    }
  }
  cout << "True: " << trueBJet << " Tagged:" << taggedBJet << " Ratio:" << (float)taggedBJet/(float)trueBJet << endl;
}

void CHiggsAnalysis::SelectGoodJets(std::vector<BaseParticle*>& goodJets, TClonesArray* branchJet, const std::vector<BaseParticle*>& goodElectrons) {
  float electronJetCloseness = 0.1;
  for (int i = 0; i < branchJet->GetEntries(); ++i) {
    BaseParticle* jet = GetJetFromBranch(i, branchJet);
    if (jet->PT() > jetMinPt && fabs(jet->Eta()) < jetMaxEta) {
      // electrons are detected as jets too
      // we want to skip these, so we will compare each jet to the previously found electrons
      // and skip it if it is one of them
      bool electronJetMatch = false;
      for (size_t ii = 0; ii < goodElectrons.size(); ++ii) {
        if (goodElectrons[ii]->P4().DeltaR(jet->P4()) < electronJetCloseness){
          electronJetMatch = true;
        }
      }

      if (!electronJetMatch) {
        goodJets.push_back(jet);
      }
    }
  }
}

void CHiggsAnalysis::ResetGoodArrays() {
}

void CHiggsAnalysis::ProcessEvent() {
  // Cutting thresholds
  float electronMinPt = 5.0, electronMaxEta = 2.5;
  float muonMinPt = 5.0, muonMaxEta = 2.5;
  std::vector<BaseParticle*> goodElectrons;
  std::vector<BaseParticle*> goodMuons;
  std::vector<BaseParticle*> goodJets6;
  std::vector<BaseParticle*> goodJets4;

  // Select good muons and electrons
  for (int i = 0; i < branchElectron->GetEntries(); ++i) {
    BaseParticle* electron = GetElectron(i);
    if (electron->PT() > electronMinPt && fabs(electron->Eta()) < electronMaxEta) {
      goodElectrons.push_back(electron);
    }
  }

  for (int i = 0; i < branchMuon->GetEntries(); ++i) {
    BaseParticle* muon = GetMuon(i);
    if (muon->PT() > muonMinPt && fabs(muon->Eta()) < muonMaxEta) {
      goodMuons.push_back(muon);
    }
  }

  Fill("LeptonMult", goodMuons.size() + goodElectrons.size());
  Fill("MuonMult", goodMuons.size());
  Fill("ElectronMult", goodElectrons.size());
  Fill("OriginalLeptonMult", branchElectron->GetEntries() + branchMuon->GetEntries());

  // Find a dileption - in this case we care about a pair of electrons or a pair of muons
  std::vector<BaseParticle*> dilepton;
  if (goodElectrons.size() == 2 && goodMuons.size() == 0) {
    dilepton.push_back(goodElectrons[0]);
    dilepton.push_back(goodElectrons[1]);
  } else if (goodElectrons.size() == 0 && goodMuons.size() == 2) {
    dilepton.push_back(goodMuons[0]);
    dilepton.push_back(goodMuons[1]);
  }

//  cout << "m: " << goodMuons.size() << " e:" << goodElectrons.size() << endl;

  Fill("Flowcut", 10);
  // Exactly two leptons
  if (dilepton.size() != 2) {
    Fill("Cut stage", 10);
    return;
  }

  Fill("Flowcut", 20);
  // Leading lepton Pt has to be > 30 GeV or leading > 20 GeV and sub-leading > 10 GeV
  if (dilepton[0]->P4().Pt() < 20.0 || (dilepton[0]->P4().Pt() < 30.0 && dilepton[1]->P4().Pt() < 10.0)) {
    Fill("Cut stage", 20);
    return;
  }

  Fill("Flowcut", 30);
  // Dileption mass has to be between 80 and 100 GeV
  float dileptonMass = (dilepton[0]->P4() + dilepton[1]->P4()).M();
  if (dileptonMass < 80.0 || dileptonMass > 100.0) {
    Fill("Cut stage", 30);
    return;
  }

  Fill("Flowcut", 40);
  // Add their mass and Pt to the histograms
  TLorentzVector tlDilepton = dilepton[0]->P4() + dilepton[1]->P4();
  // Cut Z with less then 40 GeV Pt
  if (tlDilepton.Pt() < 40.0) {
    Fill("Cut stage", 40);
    return;
  }

  Fill("Flowcut", 50);
  for (int i = 0; i < branchJet4->GetEntries(); ++i) {
    Jet* jet = (Jet*)GetJet4(i)->underlying();
    Fill("bstate4", jet->BTag);
  }

  for (int i = 0; i < branchJet6->GetEntries(); ++i) {
    Jet* jet = (Jet*)GetJet6(i)->underlying();
    Fill("bstate6", jet->BTag);
  }

  // Select good jets - different from electrons
  SelectGoodJets(goodJets4, branchJet4, goodElectrons);
  SelectGoodJets(goodJets6, branchJet6, goodElectrons);

  static int proposed = 0;
  static int tagged = 0;
  static int total = 0;
  total += goodJets6.size();
  for (size_t i = 0; i < goodJets6.size(); i++) {
    if (((Jet*)goodJets6[i])->BTag & 2) proposed++;
    if (((Jet*)goodJets6[i])->BTag & 1) tagged++;
  }

//  cout << total << "\t" << proposed << "\t" << tagged << endl;

  std::vector<BaseParticle*> boostedJets;
  //float minBoostedJetPt = 430.0;
  float minBoostedJetPt = 2.0 * 126.0 / coneSize;
  
  std::vector<BaseParticle*> resolvedJets;
  bool usedResolvedJets = false;
  if (goodJets4.size() >= 2) {
    Jet* j1 = (Jet*)goodJets4[0]->underlying();
    Jet* j2 = (Jet*)goodJets4[1]->underlying();
    if ((j1->BTag & 2) && (j2->BTag & 2)) {
      TLorentzVector bb;
      bb = goodJets4[0]->P4() + goodJets4[1]->P4();
      if (bb.M() >= 85.0 && bb.M() <= 165.0) {
	resolvedJets.push_back(goodJets4[0]);
	resolvedJets.push_back(goodJets4[1]);
	usedResolvedJets = true;
      }
    }
  } 

  if (!usedResolvedJets && goodJets6.size() > 0) {
    Jet* jet = (Jet*)goodJets6[0]->underlying();
    if (goodJets6[0]->P4().Pt() > minBoostedJetPt) {
      if (jet->BTag & 2) {
        boostedJets.push_back(goodJets6[0]);
      } else {
	Fill("Cut stage", 16);
      }
    } else {
      Fill("Cut stage", 17);
    }
  } 

  // Make histograms from what we've found
  TLorentzVector reconstructedJets;
  float btagEfficiency = 0.6;
  if (boostedJets.size() > 0) {
    reconstructedJets = boostedJets[0]->P4();
  } else if (resolvedJets.size() >= 2) {
    reconstructedJets = resolvedJets[0]->P4() + resolvedJets[1]->P4();
  } else {
    Fill("Cut stage", 60);
    return;
  }

  Fill("Flowcut", 60);
  // Filter out dijets which don't reconstruct to a standard higgs
  if (reconstructedJets.M() < 85.0 || reconstructedJets.M() > 165.0) {
    Fill("Cut stage", 70);
    return;
  }

  Fill("Flowcut", 70);
  if (boostedJets.size() > 0) {
    Fill("pT(bb) boosted", reconstructedJets.Pt(), currentXsection * btagEfficiency * btagEfficiency);
    Fill("m(bb) boosted", reconstructedJets.M(), currentXsection * btagEfficiency * btagEfficiency);
  } else if (resolvedJets.size() >= 2) {
    Fill("pT(bb) resolved", reconstructedJets.Pt(), currentXsection * btagEfficiency);
    Fill("m(bb) resolved", reconstructedJets.M(), currentXsection * btagEfficiency);
  }

  TLorentzVector higgsA = tlDilepton + reconstructedJets;
  if (useMassCut) {
    if (higgsA.M() < minMassCut || higgsA.M() > maxMassCut) {
      return;
    }
  }

  Fill("Flowcut", 80);
  if (boostedJets.size() > 0) {
    Fill("m(llbb) boosted", higgsA.M(), currentXsection * btagEfficiency * btagEfficiency);
    Fill("pT(llbb) boosted", higgsA.Pt(), currentXsection * btagEfficiency * btagEfficiency);
  } else if (resolvedJets.size() >= 2) {
    Fill("m(llbb) resolved", higgsA.M(), currentXsection * btagEfficiency);
    Fill("pT(llbb) resolved", higgsA.Pt(), currentXsection * btagEfficiency);
  }

  Fill("m(ll)", tlDilepton.M(), currentXsection);
  Fill("pT(ll)", tlDilepton.Pt(), currentXsection);

  Fill("Cut stage", -1);

  if (boostedJets.size() > 0) {
    if (higgsA.M() >200.0 && higgsA.M() < 350.0) {
      Fill("DeltaPhi(z,h)", reconstructedJets.DeltaPhi(tlDilepton));
      Fill("DeltaR(z,h)", reconstructedJets.DeltaR(tlDilepton));
      Fill("Strange Pt(z)", tlDilepton.Pt());
      Fill("Strange Pt(h)", reconstructedJets.Pt());
      Fill("Strange Eta(z)", tlDilepton.Eta());
      Fill("Strange Eta(h)", reconstructedJets.Eta());
    }
  }

  Fill("pt(ll) vs deltar(ll)", std::pair<float, float>(tlDilepton.Pt(), dilepton[0]->P4().DeltaR(dilepton[1]->P4())));
  if (resolvedJets.size() >=2 ) {
    Fill("pt(jj) vs deltar(jj)", std::pair<float, float>((resolvedJets[0]->P4()+resolvedJets[1]->P4()).Pt(), resolvedJets[0]->P4().DeltaR(resolvedJets[1]->P4())));
  }
  //ReleaseAllocatedParticles();
}

bool CHiggsAnalysis::HasHiggsParent(BaseParticle* jet) {
  float sameParticleCone = 0.3;
  for (int i = 2; i < branchParticle->GetEntries(); i++) {
    GenParticle* particle = (GenParticle*)branchParticle->At(i);
    if (jet->P4().DeltaR(particle->P4()) < sameParticleCone) {
      while (particle && particle->M1 >0) {
        if (particle->PID == 25) {
          return true;;
        }
        
        particle = (GenParticle*)branchParticle->At(particle->M1);
      }
    }
  }

  return false;
}

int CHiggsAnalysis::MaxParticleForJet(BaseParticle* jet, bool& hasBq, bool& hasHiggs) {
  int maxParticleCode = -1;
  float coneSizeForBTag = 0.3;
  hasHiggs = hasBq = false;

  for (int i = 2; i < branchParticle->GetEntries(); i++) {
    GenParticle* particle = (GenParticle*)branchParticle->At(i);
    if (particle->Status != 3) {
      continue;
    }

    if (fabs(particle->P4().Eta()) > 2.5 || particle->P4().Pt() < 1.0) {
      continue;
    }
    
    int initialParticleId = abs(particle->PID);
    if (initialParticleId == 21) {
      initialParticleId = 0;
    }

    if (jet->P4().DeltaR(particle->P4()) < coneSizeForBTag) {
      if (initialParticleId == 5) {
        hasBq = true;
      }

      if (maxParticleCode < initialParticleId) {
        maxParticleCode = initialParticleId;
      }
    }
  }

  return maxParticleCode;
}

CHiggsAnalysisNoCuts::CHiggsAnalysisNoCuts() {
  SetHistogramMinMax("1st j eta", -5, 5);
  SetHistogramMinMax("dijet Mass", 0, 2000);
}

CHiggsAnalysisNoCuts::~CHiggsAnalysisNoCuts() {
}

void CHiggsAnalysisNoCuts::ProcessEvent() {
}

void CHiggsAnalysis::ProcessHistograms() {
  CAnalysisData::ProcessHistograms();

  for (std::map<std::string, TH1F*>::iterator i = histograms.begin(); i != histograms.end(); ++i) {
    i->second->Scale(1/scaleFactor);
    i->second->Scale(20000.0);
  }
  // SetHistogramMinMax("Boosted jets mass m(H)", 40, 180);
  // SetHistogramMinMax("Resolved jets mass m(H)", 40, 180);
  // SetHistogramMinMax("Boosted jets m(A)", 500, 1100);
  // SetHistogramMinMax("Resolved jets m(A)", 500, 1100);
}

int HiggsAnalysis(const char* inputFolder, const char* inputFile, const char* outputFolder, float xsection, float coneSize) {
  CHiggsAnalysis* pha = new CHiggsAnalysis();
  pha->Initialize(inputFolder, inputFile, xsection, outputFolder, "res_with_cuts_");
  pha->SetSpecialHistogramRanges();
  pha->coneSize = coneSize;
  pha->useMassCut = false;
  std::cout << "Read " << pha->eventCount << " events." << std::endl;
  
  pha->IterateOverEvents();
  pha->SaveResults();
  
  delete pha;

  return 0;
}

int HiggsAnalysisExtraCuts(const char* inputFolder, const char* inputFile, const char* outputFolder, float xsection, float amass, float coneSize) {
  CHiggsAnalysis* pha = new CHiggsAnalysis();
  pha->Initialize(inputFolder, inputFile, xsection, outputFolder, "res_with_xtra_cuts_");
  pha->SetSpecialHistogramRanges();
  pha->coneSize = coneSize;
  pha->useMassCut = true;
  pha->minMassCut = amass - 0.1 * amass;
  pha->maxMassCut = amass + 0.05 * amass;
  std::cout << "Read " << pha->eventCount << " events." << std::endl;
  
  pha->IterateOverEvents();
  pha->SaveResults();
  
  delete pha;

  return 0;
}

//int HiggsAnalysis(const char* inputFolder, const char* inputFile, const char* outputFolder) {
//  return HiggsAnalysis(inputFolder, inputFile, outputFolder, 1.0);
//}

int HiggsTruthAnalysis(const char* inputFolder, const char* inputFile, const char* outputFolder, float heavyHiggsMass) {
  CHiggsTruth* cht = new CHiggsTruth();
  cht->Initialize(inputFolder, inputFile, 0.0, outputFolder, "res_truth_");
  cht->heavyHiggsMass = heavyHiggsMass;
  cht->IterateOverEvents();
  cht->SaveResults();
  return 0;
}

int HiggsAnalysisNoCuts(const char* inputFolder, const char* intputFile) {
  CHiggsAnalysisNoCuts* pha = new CHiggsAnalysisNoCuts();
  pha->Initialize(inputFolder, intputFile, 1.0, inputFolder, "res_no_cuts_");
  pha->IterateOverEvents();
  pha->SaveResults();
  delete pha;
  return 0;
}

int HiggsPlusBackground(const char* inputFolder, const char* inputFile, const char* inputBackgroundFolder, const char* inputBackgroundFile) {
  CHiggsAnalysisNoCuts* cha_higgs = new CHiggsAnalysisNoCuts();
  cha_higgs->Initialize(inputFolder, inputFile, 1.0, inputFolder, "res_no_cuts_");
  cha_higgs->IterateOverEvents();
  cha_higgs->SaveResults();

  CHiggsAnalysisNoCuts* cha_back = new CHiggsAnalysisNoCuts();
  cha_back->Initialize(inputBackgroundFolder, inputBackgroundFile, 1.0, inputBackgroundFolder, "res_no_cuts_");
  cha_back->IterateOverEvents();
  cha_back->SaveResults();

  CAnalysisData::DrawBoth(cha_higgs, cha_back);

  delete cha_higgs;
  delete cha_back;

  return 0;
}

int BachgroundAnalysis(const char* inputFolder, vector<TString>& inputFile, vector<float>& xsections, const char* outputFolder, float coneSize) {
  CBackGroundAnalysis* ba = new CBackGroundAnalysis();
  ba->Initialize(inputFolder, inputFile, xsections, outputFolder, "res_back_");
  std::cout << "Read " << ba->eventCount << " events." << std::endl;
  ba->SetSpecialHistogramRanges();
  ba->coneSize = coneSize;
  ba->IterateOverEvents();
  ba->SaveResults();
  delete ba;

  return 0;
}

CBackGroundAnalysis::CBackGroundAnalysis() {
  scaleFactor = 8.0 * 50000.0;
  useMassCut = false;
}

CBackGroundAnalysis::~CBackGroundAnalysis() {
}

void CBackGroundAnalysis::Initialize(const char* inputFolder, vector<TString>& inputFile, vector<float>& xsections, const char* outputFolder, const char* outputPrefix) {
  current = 0;
  inputFileList = inputFile;
  xsectionsList = xsections;
  currentInputFolder = inputFolder;

  CHiggsAnalysis::Initialize(inputFolder, inputFileList[0], xsections[0], outputFolder, outputPrefix);
}

void CBackGroundAnalysis::IterateOverEvents() {
  do
  {
    for (Long64_t event = 0; event < eventCount; ++event) {
      treeReader->ReadEntry(event);
      ProcessEvent();
    }
  } while (loadNextFile());
}

bool CBackGroundAnalysis::loadNextFile() {
  current++;
  if (current >= inputFileList.size()) {
    return false;
  }

  inputPath = inputFileList[current];
  currentXsection = xsectionsList[current];
  inputPath.Prepend(currentInputFolder);

  ReadInputFile();
  GetBranches();

  return true;
}
