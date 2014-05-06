#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TRefArray.h"
#include "TClonesArray.h"
#include "TH2D.h"
#include "TList.h"

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "Qjets/Qjets.h"
#include "Qjets/QjetsPlugin.h"

#include "JetSub.h"
#include "TreeStruct.h"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

JetSub::JetSub(const char* inputpath, const char* outputpath) {
  //intialize input data
  inputfile = new TChain("Delphes");
  inputfile->Add(inputpath);
  reader = new ExRootTreeReader(inputfile);
  b_Jet4 = reader->UseBranch("Jet4");
  b_Jet6 = reader->UseBranch("Jet6");
  b_Jet10 = reader->UseBranch("Jet10");
  b_Particle = reader->UseBranch("Particle");
  reader->UseBranch("EFlowTrack");
  reader->UseBranch("EFlowTower");
  reader->UseBranch("EFlowMuon");
  
  //intialize output data
  outputfile = new TFile(outputpath, "RECREATE");
  outputfile->cd();
  o_Jet4 = new TClonesArray("treestruct::Jet", 1);
  o_Jet6 = new TClonesArray("treestruct::Jet", 1);
  o_Jet10 = new TClonesArray("treestruct::Jet", 1);
  writer = new TTree("JetSub", "JetSub");
  writer->Branch("Jet4", &o_Jet4);
  writer->Branch("Jet6", &o_Jet6);
  writer->Branch("Jet10", &o_Jet10);
}

JetSub::~JetSub() {
  outputfile->Close();
  delete reader;
  delete inputfile;
}

void
JetSub::Run() {
  for (Int_t event_entry = 0; event_entry < reader->GetEntries(); event_entry++) {
    reader->ReadEntry(event_entry);
    
    //clear variables lists
    pass_list4.clear();
    pass_list6.clear();
    pass_list10.clear();
    qvar4.clear();
    qvar6.clear();
    qvar10.clear();
    nvar4.clear();
    nvar6.clear();
    nvar10.clear();
    
    //run Qjet and Nsubjettiness algorithms
    if (b_Jet4->GetEntries() > 0) {
      pass_list4 = Select(b_Jet4);
      qvar4.push_back(Qjet((Jet*) b_Jet4->At(pass_list4[0])));
      nvar4.push_back(Nsubjettiness((Jet*) b_Jet4->At(pass_list4[0])));
    }
    if (b_Jet6->GetEntries() > 0) {
      pass_list6 = Select(b_Jet6);
      qvar6.push_back(Qjet((Jet*) b_Jet6->At(pass_list6[0])));
      nvar6.push_back(Nsubjettiness((Jet*) b_Jet6->At(pass_list6[0])));
    }
    if (b_Jet10->GetEntries() > 0) {
      pass_list10 = Select(b_Jet10);
      qvar10.push_back(Qjet((Jet*) b_Jet10->At(pass_list10[0])));
      nvar10.push_back(Nsubjettiness((Jet*) b_Jet10->At(pass_list10[0])));
    }    
    //write single-event data into file
    Write();
  }
  
  //save TTree
  writer->Write();
}

void
JetSub::Initialize() {
  //pruning parameters
  zcut = 0.1;
  dcut_fctr = 0.5; 

  //Qjet parameters
  exp_min = 0.0;
  exp_max = 0.0;
  rigidity = 0.1;
  truncation_fctr = 0.01;
  njets = 150;
  
  //NSubjetniss parameter
  r0 = 0.3;

}

vector<int>
JetSub::Select(TClonesArray *jets) {
  vector<int> pass_list;
  
  //select leading PT jet
  int leading_entry = 0;
  Jet* leading_jet = (Jet*) jets->At(0);
  for (Int_t jet_entry = 1; jet_entry < jets->GetEntries(); jet_entry++) {
    Jet* jet = (Jet*) jets->At(jet_entry);
    if (jet->PT > leading_jet->PT) {
      if (abs(jet->Eta) < 2.0) {
	leading_entry = jet_entry;
	leading_jet = (Jet*) jets->At(jet_entry);
      }
    }
  }
  pass_list.push_back(leading_entry);
  return  pass_list;
}

void
JetSub::Write () {
  WriteJet(o_Jet4, b_Jet4, pass_list4, qvar4, nvar4);
  WriteJet(o_Jet6 ,b_Jet6, pass_list6, qvar6, nvar6);
  WriteJet(o_Jet10, b_Jet10, pass_list10, qvar10, nvar10);
  writer->Fill();
}

void
JetSub::WriteJet(TClonesArray* o_jet, TClonesArray* b_jet, vector<int>& pass_list, vector<QjetVar>& qvar, vector<NSubVar>& nvar) {
  o_jet->Clear();
  
  for (int i = 0; i < pass_list.size(); i++) {
    if (qvar.at(i).valid) {
      Jet* jet = (Jet*) b_jet->At(pass_list[i]);
      treestruct::Jet* t_jet = (treestruct::Jet*)o_jet->ConstructedAt(i);
      t_jet->PT = jet->PT;
      t_jet->Eta = jet->Eta;
      t_jet->Phi = jet->Phi;
      t_jet->Mass = jet->Mass;
      t_jet->Tau1 = nvar[i].tau1;
      t_jet->Tau2 = nvar[i].tau2;
      t_jet->Tau3 = nvar[i].tau3;
      t_jet->Mean = qvar[i].mean;
      t_jet->Std = qvar[i].std;
      t_jet->Volality = t_jet->Std/t_jet->Mean;
    }
  }
  o_jet->Compress();
}

JetSub::NSubVar
JetSub::Nsubjettiness(Jet* jet) {
  TRefArray constit_array = jet->Constituents;
  vector<fastjet::PseudoJet> constists;
  
  //extract momentum information from indicated jet
  for (Int_t i = 0; i < jet->Constituents.GetEntriesFast(); i++) {
    TObject* object = jet->Constituents.At(i);
    if (object == 0) continue;

    TLorentzVector momentum;
    momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

    if (object->IsA() == GenParticle::Class()) {
      momentum += ((GenParticle*)object)->P4();
    } else if (object->IsA() == Track::Class()) {
      momentum += ((Track*)object)->P4();
    } else if (object->IsA() == Tower::Class()) {
      momentum += ((Tower*)object)->P4();
    } else if (object->IsA() == Muon::Class()) {
      momentum += ((Muon*)object)->P4();
    } else continue;

    constists.push_back(fastjet::PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E()));
  }
  
  fastjet::JetDefinition jet_def (fastjet::kt_algorithm, 2.0*r0);
  fastjet::ClusterSequence nsub_seq(constists, jet_def);
  
  vector <double> tau(3, 0.0); 
  for (int n = 1; n <= 3; n++) {
    if (constists.size() > n) {
      // find axes
      vector<fastjet::PseudoJet> axes = nsub_seq.exclusive_jets(n);
      // calculate with closest axis
      double tau_den = 0.0;
      double* tau_num = new double[n];
      for (int i = 0; i < constists.size(); i++){
	int j_min = 0;
	double min_r = constists[i].delta_R(axes[0]);

	for (int j = 1; j < axes.size(); j++){
	  double temp_r = constists[i].delta_R(axes[j]);
	  if (temp_r < min_r) { 
	    min_r = temp_r; 
	    j_min = j;
	  }
	}
	tau_num[j_min] += constists[i].pt() * min_r;
	tau_den += constists[i].pt() * r0;
      }
      // add pieces
      for (unsigned j = 0; j < axes.size(); j++){ 
	tau[n-1] += tau_num[j]/tau_den; 
      }
    }
  }
  
  NSubVar rec(tau[0], tau[1], tau[2]);
  return rec;
}

JetSub::QjetVar
JetSub::Qjet(Jet* jet) {
  //Qjet intialization
  QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncation_fctr);
  fastjet::JetDefinition qjet_def(&qjet_plugin);
  
  TRefArray constit_array = jet->Constituents;
  vector<fastjet::PseudoJet> constists;
  
  //extract momentum information from indicated jet
  for (Int_t i = 0; i < jet->Constituents.GetEntriesFast(); i++) {
    TObject* object = jet->Constituents.At(i);
    if (object == 0) continue;

    TLorentzVector momentum;
    momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

    if (object->IsA() == GenParticle::Class()) {
      momentum += ((GenParticle*)object)->P4();
    } else if (object->IsA() == Track::Class()) {
      momentum += ((Track*)object)->P4();
    } else if (object->IsA() == Tower::Class()) {
      momentum += ((Tower*)object)->P4();
    } else if (object->IsA() == Muon::Class()) {
      momentum += ((Muon*)object)->P4();
    } else continue;

    constists.push_back(fastjet::PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E()));
  }
  
  //Now recluster indicated jet many times using Q-jet algorithm
  vector<double> masses;
  for (Int_t j = 0; j < njets; j++) {
    fastjet::ClusterSequence qjet_seq(constists, qjet_def);
    vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(qjet_seq.inclusive_jets());
    if (inclusive_jets.size() > 0)
      masses.push_back(inclusive_jets[0].m());
  }
  
  QjetVar q;
  if (masses.size()>0) {
    q.valid = true;
    //calculate mean
    for(int i = 0; i < masses.size(); i++)
      q.mean += masses[i];
    q.mean /= masses.size();

    //calculate standard deviation
    Float_t var = 0;
    for(int i = 0; i < masses.size(); i++)
      var += (masses[i]-q.mean)*(masses[i]-q.mean);
    var /= masses.size();
    q.std = sqrt(var);
  }

  return q;
}


