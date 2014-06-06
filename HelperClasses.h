#ifndef _HELPERCLASSES_H
#define _HELPERCLASSES_H

#include <vector>
#include <TLorentzVector.h>
#include <TChain.h>
#include <external/ExRootAnalysis/ExRootTreeReader.h>
#include <TH1.h>
#include <TH2.h>

using namespace std;

class GlobalData {
public:
    std::map<double, double> resolvedYield;
    std::map<double, double> boostedLowYield;
    std::map<double, double> boostedHighYield;
    std::map<double, int> resolvedCount;
    std::map<double, int> boostedLowCount;
    std::map<double, int> boostedHighCount;
    GlobalData();
    ~GlobalData();
};

class back_entry {
public:
    TString file;
    float xsection;
    
    back_entry() {
    }
    
    back_entry(const back_entry& b) {
        file = b.file;
        xsection = b.xsection;
    }
    
    back_entry& operator=(const back_entry& b) {
        file = b.file;
        xsection = b.xsection;
        
        return *this;
    }
};

class BaseParticle {
public:
    virtual Float_t PT() = 0;
    virtual Float_t Eta() = 0;
    virtual Float_t Phi() = 0;
    virtual Int_t Charge() = 0;
    virtual TLorentzVector P4() = 0;
    virtual void* underlying() = 0;
    
    virtual ~BaseParticle() {};
};

template <class GenericParticle> class Lepton : public GenericParticle, public BaseParticle {
private:
    GenericParticle* genericParticle;
public:
    Lepton(GenericParticle* p) { genericParticle = p; }
    Lepton(Lepton& anotherLepton) { genericParticle = anotherLepton.GenericParticle; }
    
    virtual Float_t PT() { return genericParticle->PT; }
    virtual Float_t Eta() { return genericParticle->Eta; }
    virtual Float_t Phi() { return genericParticle->Phi; }
    virtual Int_t Charge() { return genericParticle->Charge; }
    virtual TLorentzVector P4() { return genericParticle->P4(); }
    virtual void* underlying() { return genericParticle; }
};

class CAnalysisData {
protected:
    TString inputPath;
    TString outputPath;
    
    TFile* outputFile;
    TChain* chainDelphes;
    
    ExRootTreeReader* treeReader;
    
    TClonesArray* branchElectron;
    TClonesArray* branchMuon;
    TClonesArray* branchJet4;
    TClonesArray* branchJet6;
    TClonesArray* branchJet10;
    TClonesArray* branchParticle;
    
    std::map<std::string, std::vector<float> > histogramValues;
    std::map<std::string, std::vector<std::pair<float, float> > > histogram2Values;
    std::map<std::string, std::vector<float> > weightValues;
    std::map<std::string, TH1F*> histograms;
    std::map<std::string, TH2F*> histograms2;
    std::vector<BaseParticle*> allocatedParticles;
    std::map<std::string, float> minHistogramValues;
    std::map<std::string, float> maxHistogramValues;
public:
    CAnalysisData();
    virtual ~CAnalysisData();
    
    void Initialize(const char* inputFolder, const char* inputFile, float xsection, const char* outputFolder, const char* );
    void SimpleInitialize(const char* inputFile, const char* outputFolder, const char* prefix);
    void ReadInputFile();
    void OpenOutputFile();
    void GetBranches();
    virtual void SaveResults();
    
    void Fill(const std::string& histName, float value);
    void Fill(const std::string& histName, float value, float weight);
    void Fill(const std::string& hist2Name, std::pair<float, float> value);
    void ReleaseAllocatedParticles();
    void SetHistogramMinMax(const std::string&hist, float min, float max);
    TH1F* GetHistogram(std::string histName);
    TH2F* GetHistogram2(std::string histName);
    static bool compareFloatPairFirst(const std::pair<float, float>& p1, const std::pair<float, float>& p2) { return p1.first < p2.first; }
    static bool compareFloatPairSecond(const std::pair<float, float>& p1, const std::pair<float, float>& p2) { return p1.second < p2.second; }
    
    virtual void ProcessEvent();
    virtual void IterateOverEvents();
    virtual void ProcessHistograms();
    void Process1DHistograms();
    void Process2DHistograms();
    void FillJetData(const std::string& prefix, BaseParticle* jet);
    void FillJetData(const std::string& prefix, const TLorentzVector& jetData);
    virtual BaseParticle* GetElectron(Long64_t num);
    virtual BaseParticle* GetMuon(Long64_t num);
    virtual BaseParticle* GetJet4(Long64_t num);
    virtual BaseParticle* GetJet6(Long64_t num);
    virtual BaseParticle* GetJet10(Long64_t num);
    virtual BaseParticle* GetJetFromBranch(Long64_t num, TClonesArray* branch);
    virtual BaseParticle* GetParticle(Long64_t num);
    static void DrawBoth(CAnalysisData* a1, CAnalysisData* a2);
    
    Long64_t eventCount;
    bool normalize;
    float coneSize;
    float currentXsection;
};

class CFourVectorBranch {
public:
    CFourVectorBranch(const char* prefix);
    ~CFourVectorBranch();
    void AttachToTree(TTree* tree);
    void AddEntry(const TLorentzVector& fv);
    void Reset();
    
    vector<double> *px;
    vector<double> *py;
    vector<double> *pz;
    vector<double> *pt;
    
    string prefix;
};

class CTaus {
public:
    CTaus(const char* prefix);
    ~CTaus();
    void AttachToTree(TTree* tree);
    void AddEntry(double tau1, double tau2, double tau3);
    void Reset();
    
    vector<double> *ptau1;
    vector<double> *ptau2;
    vector<double> *ptau3;
    
    string prefix;
};

class QJetsData {
public:
    QJetsData(const char* prefix);
    ~QJetsData();
    void AttachToTree(TTree* tree);
    void AddEntry(double m, double m2, double avg, double rms, double rmss);
    void Reset();
    
    vector<double> *pm;
    vector<double> *pm2;
    vector<double> *pavg;
    vector<double> *prms;
    vector<double> *prmss;
    
    string prefix;
};

template <class t>
class SimpleScalarBranch {
public:
    SimpleScalarBranch(const char* _name) {
        name = _name;
        values = new vector<t>();
    }
    
    ~SimpleScalarBranch() {
        delete values;
    }
    
    void AttachToTree(TTree* tree) {
        tree->Branch(name.c_str(), &values);
    }
    
    void AddEntry(t value) {
        values->push_back(value);
    }
    
    void Reset() {
        values->clear();
    }
    
    vector<t> *values;
    string name;
};

extern GlobalData gd;

#endif // _HELPERCLASSES_H
