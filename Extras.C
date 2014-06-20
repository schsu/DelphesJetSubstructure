#include <cstdio>
#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TH1.h>
#include <TFile.h>
#include <TLegend.h>

#include "LocalSettings.h"
#include "AtlasStyle.h"

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

void prunedGraphs(const std::string& param, const std::string& coneType) {
    TCanvas* c = new TCanvas("pruned", "pruned", 900, 700);
    SetAtlasStyle();
    char fileName[1024];
    sprintf(fileName, "%sres-mintree_jetsub_a-zh-%dGeV.root", (const char*)outputFolder, 1500);
    TFile* f = new TFile(fileName);
    TLegend* legend = new TLegend(0.7, 0.8, 0.9, 0.7);
    
    TH1F* hist = (TH1F*)f->Get((param + ",b" + coneType).c_str());
    if (hist == NULL) {
        std::cout << "Histogram not found " << fileName << std::endl;
        return;
    }

    std::string title = param + ", " + (coneType == "h" ? "R=1.2" : "R=0.8");
    title.replace(title.find("bb"), 2, "b");
    hist->SetTitle(title.c_str());
    hist->SetMarkerStyle(21);
    hist->SetMarkerColor(kRed);
    hist->GetXaxis()->SetTitle("mA, GeV");
    hist->GetYaxis()->SetTitle("Events           ");
    legend->AddEntry(hist, "before trimming", "lp");
    hist->Draw("");
    
    TH1F* histp = (TH1F*)f->Get((param + ",b" + coneType + "p").c_str());
    if (histp == NULL) {
        std::cout << "Histogram not found " << fileName << std::endl;
        return;
    }
    
    histp->SetMarkerStyle(22);
    histp->SetMarkerColor(kBlue);
    histp->GetXaxis()->SetTitle("mA, GeV");
    histp->GetYaxis()->SetTitle("Events            ");
    legend->AddEntry(histp, "after trimming", "lp");
    histp->Draw("same");
    
    legend->Draw("ACp");
    c->Print(outputFolder + param.c_str() + coneType.c_str() + ".pdf");
    
    delete f;
    delete c;
}

int main(int argc, char* argv[]) {
    overlapmbb("m(bb),r");
    overlapmbb("m(bb),bl");
    overlapmbb("m(bb),bh");
    
    prunedGraphs("m(bb)", "l");
    prunedGraphs("m(bb)", "l");
    prunedGraphs("m(bb)", "h");
    prunedGraphs("pT(bb)", "l");
    prunedGraphs("pT(bb)", "h");
    return 0;
}
