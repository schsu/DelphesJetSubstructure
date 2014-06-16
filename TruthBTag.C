#include <iostream>
#include <memory>
#include <iomanip>

#include <TChain.h>
#include <TClonesArray.h>

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
	std::cout << "Usage: truthbtag rootfile.root" << std::endl;
	return 1;
    }

    // Reads the input root file - creates a TChain and adds the file with the 
    // name passed from the command line to it
    std::unique_ptr<TChain> delphesChain(new TChain("Delphes"));
    delphesChain->Add(argv[1]);

    // Create a tree reader and get the Particle branch
    std::unique_ptr<ExRootTreeReader> treeReader(new ExRootTreeReader(delphesChain.get()));
    // Note: ExRootTreeReader manages the memory for the branches in use. Thus, you don't have to call delete on the pointer returned by UseBranch()
    TClonesArray* branchParticle = treeReader->UseBranch("Particle");
    TClonesArray* branchJet = treeReader->UseBranch("Jet4"); // Default Delphes card generates a branch "Jet". Our analysis has three different branchs - "Jet4", "Jet6" and "Jet10"
    treeReader->UseBranch("EFlowTrack");
    treeReader->UseBranch("EFlowTower");
    treeReader->UseBranch("EFlowMuon");

    std::cout << "Using branch Particle with " << treeReader->GetEntries() << " events." << std::endl;

    size_t allJets = 0;
    size_t delphesBTaggedJets = 0;
    size_t truthBTaggedJets = 0;

    // Iterate over all events in the branch
    for (Long64_t event = 0; event < treeReader->GetEntries(); ++event) {
	// Read event number 'event'. This will update the branchParticle variable with the data for this event.
	treeReader->ReadEntry(event);

	// Iterate over all jets in the event
	for (Long64_t jnum = 0; jnum < branchJet->GetEntries(); ++jnum) {
	    // Get the jet
	    Jet* j = reinterpret_cast<Jet*>(branchJet->At(jnum));

	    double btagConeSize = 0.3;
	    bool btag = false;
	    // Now iterate over all particles for this event
	    // We can skip particles 0 and 1 - they should be the two gluons colliding
	    for (Long64_t pnum = 2; pnum < branchParticle->GetEntries(); ++pnum) {
		// Get the pnum particle from the event.
		GenParticle *p = reinterpret_cast<GenParticle*>(branchParticle->At(pnum));

		// We are interested only in particles with status=3
		if (p->Status != 3) {
		    continue;
		}

		// Skip particles with large Eta (too far away from the center of the detector)
		// or too low transverse momentum. These are adjustable in the Delphes card.
		if (j->P4().Eta() > 2.5 || j->P4().Pt() < 1.0) {
		    continue;
		}

		// If the particle is close enough to the jet and is a b/~b, we can btag the jet
		// The cone size is also adjustable in the Delphes card.
		int pdgCode = std::abs(p->PID);
		if (j->P4().DeltaR(p->P4()) < btagConeSize) {
		    if (pdgCode == 5) {
			btag = true;
			break;
		    }
		}
	    }

	    allJets++;
	    if (btag) {
		truthBTaggedJets++;
	    }

	    if (j->BTag & 1) {
		delphesBTaggedJets++;
	    }
	}
    }

    std::cout << "Total jets:           " << std::right << std::setw(10) << allJets << std::endl;
    std::cout << "B-tagged by Delphes:  " << std::right << std::setw(10) << delphesBTaggedJets << std::endl;
    std::cout << "Truth B-tagged:       " << std::right << std::setw(10) << truthBTaggedJets << std::endl;

    return 0;
}
