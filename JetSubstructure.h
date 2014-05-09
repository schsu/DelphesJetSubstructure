#ifndef _JETSUBSTRUCTURE_H_
#define _JETSUBSTRUCTURE_H_

#include <vector>
#include <TROOT.h>
#include "classes/DelphesClasses.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace std;

#include "HelperClasses.h"

class CJetSubstructure : public CAnalysisData {
private:
	double electronMinPt, electronMaxEta;
	double muonMinPt, muonMaxEta;
	double jetMinPt;
	double jetMaxEta;
	double coneSize;
	double qj_m;
	double qj_m2;
	double qj_mAverage;
	double qj_mRMSS;
	double qj_mRMS;

	CFourVectorBranch *pjets_4;
	CFourVectorBranch *pjets_6;
	CFourVectorBranch *pjets_10;
	CFourVectorBranch *pmuons;
	CFourVectorBranch *pelectrons;
	CTaus *ptaus_4;
	CTaus *ptaus_6;
	CTaus *ptaus_10;
	QJetsData *qjdata;
	QJetsData *qjfdata;
	SimpleScalarBranch<int> *btags_4;
	SimpleScalarBranch<int> *btags_6;
	SimpleScalarBranch<int> *btags_10;

	TTree* ptree;
public:
	CJetSubstructure();
	virtual ~CJetSubstructure();

	virtual void ProcessEvent();
	virtual void ProcessHistograms();
  	void SelectGoodJets(std::vector<BaseParticle*>& goodJets, TClonesArray* branchJets, const std::vector<BaseParticle*>& goodElectrons);
	void Initialize(const char* inputFolder, const char* inputFile, double xsection, double coneSize, const char* outputFolder, const char* prefix);
 	virtual void SaveResults();
 	void Reset();
 	double NSubJettiness(Jet* jet, int n, double r0);
 	void StoreJets(vector<BaseParticle*>& jets, CFourVectorBranch* pfv, CTaus* ptaus, SimpleScalarBranch<int>* btags, bool doQjets);
	template <class QJPlugin> void QJetsGeneric(Jet* jet);
 	void CalculateVolatility(Jet* jet, std::vector<TLorentzVector>& jets);
	void VerifyConstituents();
};

int JetSubstructure(const char* inputFolder, const char* inputFile, const char* outputFolder, double xsection, double coneSize);


#endif // _JETSUBSTRUCTURE_H_
