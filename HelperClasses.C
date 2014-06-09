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

using namespace std;

GlobalData gd;

void CAnalysisData::ProcessEvent() {
}

void CAnalysisData::SetHistogramMinMax(const std::string&hist, float min, float max) {
  maxHistogramValues[hist] = max;
  minHistogramValues[hist] = min;
}

void CAnalysisData::IterateOverEvents() {
  for (Long64_t event = 0; event < eventCount; ++event) {
    treeReader->ReadEntry(event);
    ProcessEvent();
  }
}

BaseParticle* CAnalysisData::GetElectron(Long64_t num) {
  BaseParticle* bl = new Lepton<Electron>((Electron*)branchElectron->At(num));
  allocatedParticles.push_back(bl);
  return bl;
}

BaseParticle* CAnalysisData::GetMuon(Long64_t num) {
  BaseParticle* bl = new Lepton<Muon>((Muon*)branchMuon->At(num));
  allocatedParticles.push_back(bl);
  return bl;
}

BaseParticle* CAnalysisData::GetJetFromBranch(Long64_t num, TClonesArray* branch) {
  BaseParticle* bl = new Lepton<Jet>((Jet*)branch->At(num));
  allocatedParticles.push_back(bl);
  return bl;
}

BaseParticle* CAnalysisData::GetJet4(Long64_t num) {
  BaseParticle* bl = new Lepton<Jet>((Jet*)branchJet4->At(num));
  allocatedParticles.push_back(bl);
  return bl;
}

BaseParticle* CAnalysisData::GetJet6(Long64_t num) {
  BaseParticle* bl = new Lepton<Jet>((Jet*)branchJet6->At(num));
  allocatedParticles.push_back(bl);
  return bl;
}

BaseParticle* CAnalysisData::GetJet10(Long64_t num) {
  BaseParticle* bl = new Lepton<Jet>((Jet*)branchJet10->At(num));
  allocatedParticles.push_back(bl);
  return bl;
}

BaseParticle* CAnalysisData::GetParticle(Long64_t num) {
  BaseParticle* bl = new Lepton<GenParticle>((GenParticle*)branchParticle->At(num));
  allocatedParticles.push_back(bl);
  return bl;
}

void CAnalysisData::ReleaseAllocatedParticles() {
  while (!allocatedParticles.empty()) {
    delete allocatedParticles.back();
    allocatedParticles.pop_back();
  }
}

CAnalysisData::CAnalysisData() {
  outputFile = NULL;
  chainDelphes = NULL;
  treeReader = NULL;
  branchMuon = NULL;
  branchElectron = NULL;

  eventCount = 0;
  normalize = false;
}

CAnalysisData::~CAnalysisData() {
  ReleaseAllocatedParticles();
  delete outputFile;
  delete chainDelphes;
  delete treeReader;
}

void CAnalysisData::SimpleInitialize(const char* inputFile, const char* outputFolder, const char* outputPrefix) {
  outputPath = inputFile;
  outputPath.Prepend(outputPrefix);
  outputPath.Prepend(outputFolder);
  OpenOutputFile();
}

void CAnalysisData::Initialize(const char* inputFolder, const char* inputFile, float xsection, const char* outputFolder, const char* outputPrefix) {
  outputPath = inputFile;
  outputPath.Prepend(outputPrefix);
  outputPath.Prepend(outputFolder);
  currentXsection = xsection;

  inputPath = inputFile;
  inputPath.Prepend(inputFolder);

  ReadInputFile();
  OpenOutputFile();
  GetBranches();
}

void CAnalysisData::ReadInputFile() {
  delete chainDelphes;
  chainDelphes = new TChain("Delphes");
  chainDelphes->Add(inputPath);
  cout << inputPath << endl;
}

void CAnalysisData::OpenOutputFile() {
  outputFile = new TFile(outputPath, "recreate");
}

void CAnalysisData::GetBranches() {
  treeReader = new ExRootTreeReader(chainDelphes);
  branchElectron = treeReader->UseBranch("Electron");
  branchMuon = treeReader->UseBranch("Muon");
  branchJet4 = treeReader->UseBranch("Jet4");
  branchJet6 = treeReader->UseBranch("Jet6");
  branchJet10 = treeReader->UseBranch("Jet10");
  branchParticle = treeReader->UseBranch("Particle");
  // Do not remove these!!
  // They are not read directly but their contents is needed to have the correct jet constituents
  treeReader->UseBranch("EFlowTrack");
  treeReader->UseBranch("EFlowTower");
  treeReader->UseBranch("EFlowMuon");

  eventCount = treeReader->GetEntries();
}

void CAnalysisData::SaveResults() {
  ProcessHistograms();
  outputFile->Write();
  outputFile->Close();
  std::cout << "Wrote " << outputPath << std::endl;
}

void CAnalysisData::Fill(const std::string& histName, float value) {
  Fill(histName, value, 1.0);
}

void CAnalysisData::Fill(const std::string& histName, float value, float weight) {
  std::map<std::string, std::vector<float> >::iterator hist = histogramValues.find(histName);
  if (hist == histogramValues.end()) {
    std::vector<float> v;
    histogramValues[histName] = v;
  }

  map<string, vector<float> >::iterator w = weightValues.find(histName);
  if (w == weightValues.end()) {
    vector<float> v;
    weightValues[histName] = v;
  }

  histogramValues[histName].push_back(value);
  weightValues[histName].push_back(weight);
}

void CAnalysisData::Fill(const std::string& hist2Name, std::pair<float, float> value) {
  std::map<std::string, std::vector<std::pair<float, float> > >::iterator hist = histogram2Values.find(hist2Name);
  if (hist == histogram2Values.end()) {
    std::vector<std::pair<float, float> > v;
    histogram2Values[hist2Name] = v;
  }
  
  histogram2Values[hist2Name].push_back(value);
}

void CAnalysisData::ProcessHistograms() {
  Process1DHistograms();
  Process2DHistograms();
}

void CAnalysisData::Process1DHistograms() {
  for (std::map<std::string, std::vector<float> >::iterator it = histogramValues.begin(); it != histogramValues.end(); ++it) {
    float minHistVal = *(std::min_element(it->second.begin(), it->second.end()));
    float maxHistVal = *(std::max_element(it->second.begin(), it->second.end()));

    bool useFixedValues = false;

    if (minHistogramValues.count(it->first) > 0) {
      minHistVal = minHistogramValues[it->first];
      maxHistVal = maxHistogramValues[it->first];
      cout << "Special max/min" << minHistVal << " " << maxHistVal << endl;
      useFixedValues = true;
    }

    float average = std::accumulate(it->second.begin(), it->second.end(), 0.0) / (it->second.size() > 0 ? it->second.size() : 1);
    std::cout << it->first.c_str() << " " << maxHistVal << " " << minHistVal << " " << average << " " << " n=" << it->second.size() << std::endl;

    float expansionFactor = 0.1;
    if (useFixedValues) {
      expansionFactor = 0;
    }

    float range = maxHistVal - minHistVal;
    TH1F* hist = new TH1F(it->first.c_str(), it->first.c_str(), 100, minHistVal - expansionFactor * range, maxHistVal + expansionFactor * range);
    vector<float> weights = weightValues[it->first];
    if (weights[0] != 1.0) {
      hist->Sumw2();
    }
  
    for (size_t i = 0; i < it->second.size(); ++i) {
      hist->Fill(it->second[i], weights[i]);
    }

    histograms[it->first] = hist;
  }
}

void CAnalysisData::Process2DHistograms() {
  for (std::map<std::string, std::vector<std::pair<float, float> > >::iterator it = histogram2Values.begin(); it != histogram2Values.end(); ++it) {
    float minFirstVal = std::min_element(it->second.begin(), it->second.end(), compareFloatPairFirst)->first;
    float maxFirstVal = std::max_element(it->second.begin(), it->second.end(), compareFloatPairFirst)->first;
    float minSecondVal = std::min_element(it->second.begin(), it->second.end(), compareFloatPairSecond)->second;
    float maxSecondVal = std::max_element(it->second.begin(), it->second.end(), compareFloatPairSecond)->second;
    
    float expansionFactor = 0.1;
    float rangex = maxFirstVal - minFirstVal;
    float rangey = maxSecondVal - minSecondVal;
    TH2F* hist = new TH2F(it->first.c_str(), it->first.c_str(), 
			  100, minFirstVal - expansionFactor * rangex, maxFirstVal + expansionFactor * rangex, 
			  100, minSecondVal - expansionFactor * rangey, maxSecondVal + expansionFactor * rangey);
    for (size_t i = 0; i < it->second.size(); ++i) {
      hist->Fill(it->second[i].first, it->second[i].second);
    }

    histograms2[it->first] = hist;
  }
}

void CAnalysisData::DrawBoth(CAnalysisData* d1, CAnalysisData* d2) {
  for (std::map<std::string, TH1F*>::iterator it = d1->histograms.begin(); it != d1->histograms.end(); it++) {
    std::map<std::string, TH1F*>::iterator other = d2->histograms.find(it->first);
    if (other != d2->histograms.end()) {
      TCanvas* c = new TCanvas(it->first.c_str(), "h", 800, 600);

      Float_t maxval = it->second->GetMaximum() * 1.1;
      Float_t scale = gPad->GetUymax()/maxval;
      cout << scale << endl;
      it->second->Scale(scale);
      it->second->Draw();

      maxval = other->second->GetMaximum() * 1.1;
      scale = gPad->GetUymax()/maxval;
      other->second->SetLineColor(kRed);
      cout << scale << endl;
      other->second->Scale(scale);
      other->second->Draw("same");
      c->Update();

      TImage* img = TImage::Create();
      img->FromPad(c);
      TString imageName = it->first.c_str();
      imageName.Append(".png");
      img->WriteImage(imageName);
      delete img;
      delete c;
    }
  }
}

void CAnalysisData::FillJetData(const std::string& prefix, BaseParticle* jet) {
  FillJetData(prefix, jet->P4());
}

void CAnalysisData::FillJetData(const std::string& prefix, const TLorentzVector& jetData) {
  Fill((prefix + "Pt").c_str(), jetData.Pt());
  Fill((prefix + "Eta").c_str(), jetData.Eta());
  Fill((prefix + "Mass").c_str(), jetData.M());
}

TH1F* CAnalysisData::GetHistogram(std::string histName) {
  std::map<std::string, TH1F*>::iterator hist = histograms.find(histName);
  if (hist == histograms.end()) {
    return NULL;
  }

  return hist->second;
}

TH2F* CAnalysisData::GetHistogram2(std::string histName) {
  std::map<std::string, TH2F*>::iterator hist = histograms2.find(histName);
  if (hist == histograms2.end()) {
    return NULL;
  }

  return hist->second;
}

CFourVectorBranch::CFourVectorBranch(const char* pref) {
  prefix = pref;
  px = new vector<double>();
  py = new vector<double>();
  pz = new vector<double>();
  pt = new vector<double>();
}

CFourVectorBranch::~CFourVectorBranch() {
  delete px;
  delete py;
  delete pz;
  delete pt;
}

void
CFourVectorBranch::AttachToTree(TTree* tree) {
  tree->Branch((prefix + "_x").c_str(), &px);
  tree->Branch((prefix + "_y").c_str(), &py);
  tree->Branch((prefix + "_z").c_str(), &pz);
  tree->Branch((prefix + "_t").c_str(), &pt);
}

void 
CFourVectorBranch::AddEntry(const TLorentzVector& fv) {
  px->push_back(fv.X());
  py->push_back(fv.Y());
  pz->push_back(fv.Z());
  pt->push_back(fv.T());
}

void
CFourVectorBranch::Reset() {
  px->clear();
  py->clear();
  pz->clear();
  pt->clear();
}

CTaus::CTaus(const char* pref) {
  prefix = pref;
  ptau1 = new vector<double>();
  ptau2 = new vector<double>();
  ptau3 = new vector<double>();
}

CTaus::~CTaus() {
  delete ptau1;
  delete ptau2;
  delete ptau3;
}

void
CTaus::AttachToTree(TTree* tree) {
  tree->Branch(("tau1_" + prefix).c_str(), &ptau1);
  tree->Branch(("tau2_" + prefix).c_str(), &ptau2);
  tree->Branch(("tau3_" + prefix).c_str(), &ptau3);
}

void
CTaus::AddEntry(double tau1, double tau2, double tau3) {
  ptau1->push_back(tau1);
  ptau2->push_back(tau2);
  ptau3->push_back(tau3);
}

void
CTaus::Reset() {
  ptau1->clear();
  ptau2->clear();
  ptau3->clear();
}

QJetsData::QJetsData(const char* pref) {
  prefix = pref;
  pm = new vector<double>();
  pm2 = new vector<double>();
  pavg = new vector<double>();
  prms = new vector<double>();
  prmss = new vector<double>();
}

QJetsData::~QJetsData() {
  delete pm;
  delete pm2;
  delete pavg;
  delete prmss;
  delete prms;
}

void
QJetsData::AttachToTree(TTree* tree) {
  tree->Branch((prefix+"qj_m").c_str(), &pm);
  tree->Branch((prefix+"qj_m2").c_str(), &pm2);
  tree->Branch((prefix+"qj_avg").c_str(), &pavg);
  tree->Branch((prefix+"qj_rms").c_str(), &prms);
  tree->Branch((prefix+"qj_rmss").c_str(), &prmss);
}

void
QJetsData::AddEntry(double m, double m2, double avg, double rms, double rmss) {
  pm->push_back(m);
  pm2->push_back(m2);
  pavg->push_back(avg);
  prms->push_back(rms);
  prmss->push_back(rmss);
}

void
QJetsData::Reset() {
  pm->clear();
  pm2->clear();
  pavg->clear();
  prms->clear();
  prmss->clear();
}

GlobalData::GlobalData() {
}

GlobalData::~GlobalData() {
    if (resolvedYield.size() > 0) {

	std::cout << "Resolved yields:" << std::endl;
	for (const std::pair<double, double>& p : resolvedYield) {
	    std::cout << p.first << "->" << p.second << std::endl;
	}

	std::cout << "Boosted low yields:" << std::endl;
	for (const std::pair<double, double>& p : boostedLowYield) {
	    std::cout << p.first << "->" << p.second << std::endl;
	}

	std::cout << "Boosted high yields: " << std::endl;
	for (const std::pair<double, double>& p : boostedHighYield) {
	    std::cout << p.first << "->" << p.second << std::endl;
	}
    }
}
