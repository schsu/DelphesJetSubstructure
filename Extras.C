#include <cstdio>
#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TH1.h>
#include <TFile.h>

#include "LocalSettings.h"

void overlapmbb(const std::string& channel) {
    TCanvas* c = new TCanvas("overlapmbb", "overlapmbb", 900, 700);
    for (int i = 0; i < amassesCount; ++i) {
        char fileName[1024];
        std::sprintf(fileName, "%sres-mintree_jetsub_a-zh-%dGeV.root", (const char*)outputFolder, amasses[i]);
        std::cout << fileName << std::endl;
        TFile* f = new TFile(fileName);
        TH1F* hist = (TH1F*)f->Get(channel.c_str());
        if (hist == NULL) {
            std::cout << "Histogram not found " << fileName << std::endl;
            continue;
        }
        
        hist->Draw(i > 0 ? "same" : "");
    }

    c->Print(outputFolder + "overlapmbb" + channel.c_str() + ".pdf");
    delete c;
}

int main(int argc, char* argv[]) {
    overlapmbb("m(bb),r");
    overlapmbb("m(bb),bl");
    overlapmbb("m(bb),bh");
    return 0;
}
