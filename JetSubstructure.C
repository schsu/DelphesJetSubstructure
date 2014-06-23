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

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"

#include "qjet_old/QjetsPlugin.h"
#include "qjet_old/Qjets.C"
#include "qjet_old/QjetsPlugin.C"

#include "QjetsPlugin.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <numeric>
#include <fstream>
#include <chrono>

#include "HelperClasses.h"
#include "JetSubstructure.h"

int counter = 0;

CJetSubstructure::CJetSubstructure() {
    electronMinPt = 5.0;
    electronMaxEta = 2.5;
    muonMinPt = 5.0;
    muonMaxEta = 2.5;
    jetMinPt = 20.0;
    jetMaxEta = 2.5;
    
    ptree = NULL;
    pmuons = NULL;
    pelectrons = NULL;
    pjets_4 = NULL;
    pjets_6 = NULL;
    pjets_10 = NULL;
    prjets_6 = NULL;
    prjets_10 = NULL;
    ptaus_4 = NULL;
    ptaus_6 = NULL;
    ptaus_10 = NULL;
    qjdata = NULL;
    qjfdata = NULL;
    btags_4 = NULL;
    btags_6 = NULL;
    btags_10 = NULL;
    
    counter = 0;
}

CJetSubstructure::~CJetSubstructure() {
}

void
CJetSubstructure::ProcessEvent() {
    counter++;
    if ((counter % 1000) == 0) {
        std::cout << "#";
        std::cout.flush();
    }
    
    std::vector<BaseParticle*> goodElectrons;
    std::vector<BaseParticle*> goodMuons;
    std::vector<BaseParticle*> goodJets6;
    std::vector<BaseParticle*> goodJets4;
    std::vector<BaseParticle*> goodJets10;
    
    Reset();
	
    // Select good muons and electrons (meeting the minPt/maxEta restrictions)
    for (int i = 0; i < branchElectron->GetEntries(); ++i) {
        BaseParticle* electron = GetElectron(i);
        pelectrons->AddEntry(electron->P4());
        if (electron->PT() > electronMinPt && fabs(electron->Eta()) < electronMaxEta) {
            goodElectrons.push_back(electron);
        }
    }
    
    for (int i = 0; i < branchMuon->GetEntries(); ++i) {
        BaseParticle* muon = GetMuon(i);
        pmuons->AddEntry(muon->P4());
        if (muon->PT() > muonMinPt && fabs(muon->Eta()) < muonMaxEta) {
            goodMuons.push_back(muon);
        }
    }
    
    // Select good jets - different from electrons
    SelectGoodJets(goodJets4, branchJet4, goodElectrons);
    SelectGoodJets(goodJets6, branchJet6, goodElectrons);
    SelectGoodJets(goodJets10, branchJet10, goodElectrons);
    
    StoreJets(goodJets4, pjets_4, NULL, ptaus_4, btags_4, 0.4, false);
    StoreJets(goodJets6, pjets_6, prjets_6, ptaus_6, btags_6, 0.8, false);
    StoreJets(goodJets10, pjets_10, prjets_10, ptaus_10, btags_10, 1.2, false);
    
    ptree->Fill();
}

void
CJetSubstructure::StoreJets(vector<BaseParticle*>& jets, CFourVectorBranch* pfv, CFourVectorBranch* prfv, CTaus* ptaus, SimpleScalarBranch<int>* btags, double cone, bool doQjets) {
	for (size_t i = 0; i < jets.size(); ++i) {
        Jet* ji = (Jet*)jets[i]->underlying();
        pfv->AddEntry(jets[i]->P4());
        double tau1 = NSubJettiness(ji, 1, 0.3);
        double tau2 = NSubJettiness(ji, 2, 0.3);
        double tau3 = NSubJettiness(ji, 3, 0.3);
        ptaus->AddEntry(tau1, tau2, tau3);
        btags->AddEntry(ji->BTag);
        
        // Do some prunning
        if (prfv != NULL) {
	    TLorentzVector prunedJet = Prune(ji, cone);
            prfv->AddEntry(prunedJet);
        }
        
        if (!doQjets || i > 0) {
            qjdata->AddEntry(0.0, 0.0, 0.0, 0.0, 0.0);
            qjfdata->AddEntry(0.0, 0.0, 0.0, 0.0, 0.0);
        } else {
            std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
            QJetsGeneric<QjetsPlugin>((Jet*)jets[0]->underlying());
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> totalTime = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
            double oldTime = totalTime.count();
            Fill("OldQJets (s)", totalTime.count());
            qjdata->AddEntry(qj_m, qj_m2, qj_mAverage, qj_mRMS, qj_mRMSS);
            
            start = std::chrono::high_resolution_clock::now();
            QJetsGeneric<fastqjets::QjetsPlugin>((Jet*)jets[0]->underlying());
            end = std::chrono::high_resolution_clock::now();
            totalTime = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
            Fill("NewQJets (s)", totalTime.count());
            double newTime = totalTime.count();
            qjfdata->AddEntry(qj_m, qj_m2, qj_mAverage, qj_mRMS, qj_mRMSS);
            if (oldTime > 0.0 && newTime > 0.0) {
                Fill("Ratio (s/s)", newTime/oldTime);
            }
        }
	}
}

void
CJetSubstructure::SelectGoodJets(std::vector<BaseParticle*>& goodJets, TClonesArray* branchJet, const std::vector<BaseParticle*>& goodElectrons) {
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

void
CJetSubstructure::Initialize(const char* inputFolder, const char* inputFile, double xsection, double coneSizeParam, const char* outputFolder, const char* prefix) {
	CAnalysisData::Initialize(inputFolder, inputFile, xsection, outputFolder, prefix);
	coneSize = coneSizeParam;
    
	ptree = new TTree("DelphesNTup", "DelphesNTup");
	// Create the branches
	pjets_4 = new CFourVectorBranch("jets_antikt_4");
	pjets_4->AttachToTree(ptree);
    
	pjets_6 = new CFourVectorBranch("jets_antikt_6");
	pjets_6->AttachToTree(ptree);
    
	pjets_10 = new CFourVectorBranch("jets_antikt_10");
	pjets_10->AttachToTree(ptree);
    
    prjets_6 = new CFourVectorBranch("jets_pruned_6");
    prjets_6->AttachToTree(ptree);
    
    prjets_10 = new CFourVectorBranch("jets_pruned_10");
    prjets_10->AttachToTree(ptree);
    
	ptaus_4 = new CTaus("4");
	ptaus_4->AttachToTree(ptree);
    
	ptaus_6 = new CTaus("6");
	ptaus_6->AttachToTree(ptree);
    
	ptaus_10 = new CTaus("10");
	ptaus_10->AttachToTree(ptree);
    
	qjdata = new QJetsData("old");
	qjdata->AttachToTree(ptree);
    
	qjfdata = new QJetsData("new");
	qjfdata->AttachToTree(ptree);
    
	btags_4 = new SimpleScalarBranch<int>("btags4");
	btags_4->AttachToTree(ptree);
    
	btags_6 = new SimpleScalarBranch<int>("btags6");
	btags_6->AttachToTree(ptree);
    
	btags_10 = new SimpleScalarBranch<int>("btags10");
	btags_10->AttachToTree(ptree);
	
	pmuons = new CFourVectorBranch("muons");
	pmuons->AttachToTree(ptree);
    
	pelectrons = new CFourVectorBranch("electrons");
	pelectrons->AttachToTree(ptree);
}

void
CJetSubstructure::SaveResults() {
	CAnalysisData::SaveResults();
}

void
CJetSubstructure::Reset() {
	pjets_4->Reset();
	pjets_6->Reset();
	pjets_10->Reset();
    prjets_6->Reset();
    prjets_10->Reset();
	pelectrons->Reset();
	pmuons->Reset();
	ptaus_4->Reset();
	ptaus_6->Reset();
	ptaus_10->Reset();
	qjdata->Reset();
	qjfdata->Reset();
	btags_4->Reset();
	btags_6->Reset();
	btags_10->Reset();
}

double
CJetSubstructure::NSubJettiness(Jet* jet, int n, double r0) {
	double tau = 0.0;
	TRefArray constituents = jet->Constituents;
    
	vector<fastjet::PseudoJet> constJets;
	for (Int_t i = 0; i < constituents.GetEntriesFast(); ++i) {
		TLorentzVector momentum;
		momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
       	TObject* object = jet->Constituents.At(i);
        
		// Check if the constituent is accessible
        if(object == 0) continue;
        
        if(object->IsA() == GenParticle::Class()) {
         	momentum += ((GenParticle*) object)->P4();
        }
        else if(object->IsA() == Track::Class()) {
          	momentum += ((Track*) object)->P4();
        }
        else if(object->IsA() == Tower::Class()) {
          	momentum += ((Tower*) object)->P4();
		}
        else if(object->IsA() == Muon::Class()) {
          	momentum += ((Muon*) object)->P4();
        }
        else {
        	continue;
        }
        
		constJets.push_back(fastjet::PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E()));
	}
    
	if (constJets.size() > n) {
		// find axes
		fastjet::JetDefinition jet_def = fastjet::JetDefinition(fastjet::kt_algorithm,2.0*r0);
		fastjet::ClusterSequence* cs = new fastjet::ClusterSequence(constJets,jet_def);
		vector<fastjet::PseudoJet> axes = cs->exclusive_jets(n);
		cs->delete_self_when_unused();
		// calculate with closest axis
		double tau_den = 0.0;
		double* tau_num = new double[n];
		for (unsigned i=0; i<constJets.size(); i++){
			int j_min=0;
			double min_r = constJets[i].delta_R(axes[0]);
            
			for (unsigned j=1; j<axes.size(); j++){
			  	double temp_r = constJets[i].delta_R(axes[j]);
			  	if (temp_r < min_r) {
			  		min_r = temp_r;
			  		j_min = j;
			  	}
			}
            
			tau_num[j_min] += constJets[i].pt() * min_r;
			tau_den += constJets[i].pt() * r0;
		}
        
		// add pieces
		for (unsigned j=0; j<axes.size(); j++){
			tau += tau_num[j]/tau_den;
		}
	}
    
	return tau;
}

void
CJetSubstructure::VerifyConstituents() {
    for (int i = 0; i < branchJet10->GetEntries(); ++i) {
        BaseParticle *bp = GetJetFromBranch(i, branchJet10);
        Jet *jet = (Jet*)bp->underlying();
        TRefArray constituents = jet->Constituents;
        fstream fs("short.bin", ios::out | ios::app | ios::binary);
        double fv[4];
        fv[0] = jet->P4().X();
        fv[1] = jet->P4().Y();
        fv[2] = jet->P4().Z();
        fv[3] = jet->P4().T();
        fs.write((const char*)&fv, sizeof(fv));
        int s = constituents.GetEntriesFast();
        fs.write((const char*)&s, sizeof(s));
        for (Int_t i = 0; i < constituents.GetEntriesFast(); ++i) {
            TLorentzVector momentum;
            momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
            TObject* object = jet->Constituents.At(i);
            if (object != 0) {
                if(object->IsA() == GenParticle::Class()) {
                    momentum += ((GenParticle*) object)->P4();
                }
                else if(object->IsA() == Track::Class()) {
                    momentum += ((Track*) object)->P4();
                }
                else if(object->IsA() == Tower::Class()) {
                    momentum += ((Tower*) object)->P4();
                }
                else if(object->IsA() == Muon::Class()) {
                    momentum += ((Muon*) object)->P4();
                }
                else {
                    //continue;
                }
            } else {
            }
            
            fv[0] = momentum.X();
            fv[1] = momentum.Y();
            fv[2] = momentum.Z();
            fv[3] = momentum.T();
            fs.write((const char*)&fv, sizeof(fv));
        }
    }
}

TLorentzVector
CJetSubstructure::Prune(Jet* jet, double cone) {
    TLorentzVector result;
    result.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    // Gather the constituents
    TRefArray constituents = jet->Constituents;
    vector<fastjet::PseudoJet> constJets;
    for (Int_t i = 0; i < constituents.GetEntriesFast(); ++i) {
        TLorentzVector momentum;
        momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        TObject* object = jet->Constituents.At(i);
        
        // Check if the constituent is accessible
        if(object == 0) continue;
        
        if(object->IsA() == GenParticle::Class()) {
            momentum += ((GenParticle*) object)->P4();
        }
        else if(object->IsA() == Track::Class()) {
            momentum += ((Track*) object)->P4();
        }
        else if(object->IsA() == Tower::Class()) {
            momentum += ((Tower*) object)->P4();
        }
        else if(object->IsA() == Muon::Class()) {
            momentum += ((Muon*) object)->P4();
        }
        else {
            continue;
        }
        
        constJets.push_back(fastjet::PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E()));
    }
	
	if (constJets.size() > 0) {
		// find axes
		fastjet::JetDefinition jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm,cone);
		fastjet::ClusterSequence cs(constJets,jet_def);
		//std::cout << cs.inclusive_jets().size() << std::endl;
        std::vector<fastjet::PseudoJet> inclusiveJets = sorted_by_pt(cs.inclusive_jets());
        if (inclusiveJets.size() > 0) {
	  //std::cout << cs.inclusive_jets()[0].pt() << " " << cs.inclusive_jets()[0].m() << std::endl;
		//std::cout << jet->P4().Pt() << " " << jet->P4().M() << std::endl << std::endl;
            fastjet::PseudoJet leadingJet = inclusiveJets[0];
            fastjet::JetAlgorithm pruneAlgo = fastjet::cambridge_algorithm;
            double zcut = 0.1;
            double rcut = 0.5;
            fastjet::Pruner pruner(pruneAlgo, zcut, rcut);
	    //std::cout << pruner.description() << std::endl;
	    fastjet::PseudoJet prunedJet = pruner(leadingJet);
            if (!prunedJet.has_structure_of<fastjet::Pruner>()) {
                std::cout << "Error: Does not have pruned jet structure!" << std::endl;
            }

	    fastjet::JetDefinition jet_def_trim = fastjet::JetDefinition(fastjet::antikt_algorithm, 0.3);
	    fastjet::Filter filter(jet_def_trim, fastjet::SelectorPtFractionMin(0.05));
	    prunedJet = filter(leadingJet);

            result.SetPxPyPzE(prunedJet.px(), prunedJet.py(), prunedJet.pz(), prunedJet.e());
        }
	}
    
    return result;
}

template <class QJPlugin> void
CJetSubstructure::QJetsGeneric(Jet* jet) {
    double zcut = 0.1f;
    double dcut = 0.5f;
    double exp_min = 0.01f;
    double exp_max = 0.01f;
    double rigidity = 0.1f;
    double njets = 75.0;
    
    QJPlugin qjPlugin(zcut, dcut, exp_min, exp_max, rigidity);
	
    fastjet::JetDefinition qjet_def(&qjPlugin);
    
    TRefArray constituents = jet->Constituents;
    vector<fastjet::PseudoJet> constJets;
    for (Int_t i = 0; i < constituents.GetEntriesFast(); ++i) {
        TLorentzVector momentum;
        momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
        TObject* object = jet->Constituents.At(i);
        
        // Check if the constituent is accessible
        if(object == 0) continue;
        
        if(object->IsA() == GenParticle::Class()) {
            momentum += ((GenParticle*) object)->P4();
        }
        else if(object->IsA() == Track::Class()) {
            momentum += ((Track*) object)->P4();
        }
        else if(object->IsA() == Tower::Class()) {
            momentum += ((Tower*) object)->P4();
        }
        else if(object->IsA() == Muon::Class()) {
            momentum += ((Muon*) object)->P4();
        }
        else {
            continue;
        }
        
        constJets.push_back(fastjet::PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E()));
    }
	
    std::vector<TLorentzVector> qjets;
    for (int i = 0; i < njets; ++i) {
        fastjet::ClusterSequence qjetSeq(constJets, qjet_def);
		
        vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjetSeq.inclusive_jets());
        if (inclusive_jets2.size() < 1)
            return;
		
        fastjet::PseudoJet QJetP = inclusive_jets2[0];
        
        TLorentzVector QJet;
        QJet.SetPtEtaPhiE(QJetP.pt(), QJetP.eta(), QJetP.phi(), QJetP.e());
        qjets.push_back(QJet);
    }
    
    CalculateVolatility(jet, qjets);
}

void
CJetSubstructure::CalculateVolatility(Jet* jet, std::vector<TLorentzVector>& jets) {
	qj_m = 0.0;
	qj_m2 = 0.0;
    
	for(size_t iJ = 0; iJ < jets.size(); iJ++){
		double mJ = jets[iJ].M();
		qj_m  += mJ;
		qj_m2 += (mJ * mJ);
	}
    
	double dNJets =  (double) jets.size();
    
	qj_mAverage = qj_m / dNJets;
	qj_mRMSS    = qj_m2 / dNJets - qj_mAverage*qj_mAverage;
    
	if(qj_mRMSS < 0.0) {
		qj_mRMS = 0.0;
	} else {
		qj_mRMS     = TMath::Sqrt(qj_mRMSS);
	}
}

int
JetSubstructure(const char* inputFolder, const char* inputFile, const char* outputFolder, double xsection, double coneSize) {
	CJetSubstructure js;
	js.Initialize(inputFolder, inputFile, xsection, coneSize, outputFolder, "mintree_jetsub_");
	js.IterateOverEvents();
	js.SaveResults();
    
	return 0;
}

void
CJetSubstructure::ProcessHistograms() {
    SetHistogramMinMax("NewQJets (s)", 0.0, 2.0);
    SetHistogramMinMax("OldQJets (s)", 0.0, 2.0);
    CAnalysisData::ProcessHistograms();
}
