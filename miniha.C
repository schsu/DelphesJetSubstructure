#include <string>
#include <iostream>
#include <vector>
#include <map>

#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TGraph.h>

#include "DelphesNTuple.h"
#include "LocalSettings.h"
#include "HiggsHist.h"
#include "HelperClasses.h"
#include "AtlasStyle.h"

std::map<double, double> rbLow;
std::map<double, double> rbHigh;
std::map<double, double> signResolved;
std::map<double, double> signHigh;
std::map<double, double> signLow;
std::map<double, double> effResolved;
std::map<double, double> effHigh;
std::map<double, double> effLow;

void
AnalyzeBackground() {
    std::vector<const char*> filenames;
    std::vector<double> xsections;
    for (int j = 0; j < 1; ++j) {
        for (int i = 1; i <= 8; ++i) {
            char buffer[2048];
            sprintf(buffer, "mintree_jetsub_background-%d.root", 100 + j*9 + i);
            filenames.push_back(strdup(buffer));
        }
        
        xsections.push_back(14.59);
        xsections.push_back(0.7155);
        xsections.push_back(0.0668);
        xsections.push_back(0.01246);
        xsections.push_back(0.003444);
        xsections.push_back(0.001215);
        xsections.push_back(0.0007277);
        xsections.push_back(0.0001713);
    }
    
    HiggsHist(inputFolder, filenames, xsections, outputFolder, 1, 50000);
    std::cout << "Done with background" << std::endl;
}

void
AnalyzeSignalForAllMasses() {
    int masses[] = { 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000 };
    double xsections[] = { 0.405, 0.493, 0.0473, 0.0114, 0.00431, 0.00205, 0.00106, 0.000577, 0.000326, 0.000189 };
    for (int i = 0; i < sizeof(masses)/sizeof(masses[0]); ++i) {
        char fileName[1024];
        sprintf(fileName, "mintree_jetsub_a-zh-%dGeV.root", masses[i]);
        HiggsHist(inputFolder, fileName, outputFolder, xsections[i], masses[i], 50000);
        std::cout << "Done with " << masses[i] << "GeV" << std::endl;
        std::cout.flush();
    }
}

void
AnalyzeBackgroundSlices() {
    double xsections[] = { 14.59, 0.7155, 0.0668, 0.01246, 0.003444, 0.001215, 0.0007277, 0.0001713 };
    double masses[] = { 0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0 };
    for (int i = 0; i < 8; ++i) {
        char fileName[1024];
        sprintf(fileName, "mintree_jetsub_background-%d.root", 101+i);
        HiggsHist(inputFolder, fileName, outputFolder, xsections[i], 1.0, 50000);
        std::cout << "Done with background " << i << std::endl;
        std::cout.flush();
    }
}

std::map<double, double>
CalculateSignificance(std::map<double, double>& yield) {
    std::map<double, double> significance;
    for (const std::pair<double, double>& p : yield) {
        if (p.first < 100) {
            continue;
        }
        
        significance[p.first] = p.second / std::sqrt(yield[1.0]);
    }
    
    return significance;
}

std::map<double, double>
CalculateEfficiency(std::map<double, int>& count) {
    std::map<double, double> efficiency;
    for (const std::pair<double, int>& p : count) {
        if (p.first < 100) {
            continue;
        }
        
        efficiency[p.first] = p.second / 50000.0f;
    }
    
    return efficiency;
}

int pass = 0;

void
PlotYields(const std::string& graphName, const std::string& name, std::map<double, double> yields, TLegend* legend, const std::string& yaxisTitle) {
    int colors[] = {kBlue, kRed, kGreen, kOrange, kBlack};
    int lastColor = sizeof(colors)/sizeof(colors[0])-1;
    double miny = std::min_element(yields.begin(), yields.end(), [](std::pair<double, double> p, std::pair<double, double> q) { return p.second < q.second; })->second;
    double maxy = std::max_element(yields.begin(), yields.end(), [](std::pair<double, double> p, std::pair<double, double> q) { return p.second < q.second; })->second;
    
    TH1F* hist = new TH1F((graphName+name).c_str(), graphName.c_str(), 1800, 200, 1100);
    hist->SetMarkerStyle(21+pass);
    hist->SetMarkerColor(pass > lastColor ? lastColor : colors[pass]);
    hist->GetXaxis()->SetTitle("mA, GeV");
    hist->GetXaxis()->SetTitleOffset(1.0);
    hist->GetYaxis()->SetTitle(yaxisTitle.c_str());
    hist->GetYaxis()->SetTitleOffset(1.0);

    for (const std::pair<double, double>& p : yields) {
        // don't plot background yields
        if (p.first < 100) {
            continue;
        }
        
        std::cout << p.first << " " << p.second << std::endl;
        hist->Fill(p.first, p.second);
    }
    
    if (pass == 0) {
        hist->Draw("p");
    } else {
        hist->Draw("samep");
    }
    
    legend->AddEntry(hist, name.c_str(), "lp");
    
    //c->Print(outputFolder + "yields.pdf");
    ++pass;
}

void
PlotSignificance() {
    pass = 0;
    TFile* signFile = new TFile(outputFolder + "sign.root", "recreate");
    TCanvas* c = new TCanvas("Significance", "Significance", 900, 600);
    c->SetLogy();
    c->Print(outputFolder + "sign.pdf[");
    TLegend* legend = new TLegend(0.7, 0.8, 0.9, 0.7);
    SetAtlasStyle();
    
    PlotYields("Significance", "resolved", signResolved, legend, "");
    PlotYields("Significance", "resolved+boosted 0.8", signLow, legend, "");
    PlotYields("Significance", "resolved+boosted 1.2", signHigh, legend, "");
    
    legend->Draw("ACp");
    c->Print(outputFolder + "sign.pdf");
    
    TCanvas* c2 = new TCanvas("Significance", "Significance", 900, 600);
    c2->Print(outputFolder + "sign.pdf]");
    signFile->Write();
    signFile->Close();
}

void
PlotYields() {
    pass = 0;
    TFile* yieldsFile = new TFile(outputFolder + "yields.root", "recreate");
    // Start the pdf file
    TCanvas * c = new TCanvas("Yields", "Yields", 900, 600);
    c->SetLogy();
    c->Print(outputFolder + "yields.pdf[");
    TLegend* legend = new TLegend(0.7, 0.8, 0.9, 0.7);
    SetAtlasStyle();
    
    PlotYields("Yields", "resolved", gd.resolvedYield, legend, "Events/GeV");
    PlotYields("Yields", "resolved+boosted 0.8", rbLow, legend, "Events/GeV");
    PlotYields("Yields", "resolved+boosted 1.2", rbHigh, legend, "Events/GeV");
    
    legend->Draw("ACp");
    c->Print(outputFolder + "yields.pdf");
    
    // Finish the pdf file
    TCanvas * c2 = new TCanvas("Yields", "Yields", 900, 00);
    c2->Print(outputFolder + "yields.pdf]");
    yieldsFile->Write();
    yieldsFile->Close();
}

void
PlotEfficiency() {
    pass = 0;
    TFile* efficiencyFile = new TFile(outputFolder + "efficiency.root", "recreate");
    TCanvas * c = new TCanvas("Efficiency", "Efficiency", 900, 600);
    c->Print(outputFolder + "efficiency.pdf[");
    TLegend* legend = new TLegend(0.2, 0.8, 0.4, 0.7);
    SetAtlasStyle();
    
    PlotYields("Efficiency", "resolved", effResolved, legend, "");
    PlotYields("Efficiency", "boosted 1.2", effHigh, legend, "");
    PlotYields("Efficiency", "boosted 0.8", effLow, legend, "");
    
    legend->Draw("ACp");
    c->Print(outputFolder + "efficiency.pdf");
    
    TCanvas* c2 = new TCanvas("Efficiency", "Efficiency", 100, 100);
    c2->Print(outputFolder + "efficiency.pdf]");
    efficiencyFile->Write();
    efficiencyFile->Close();
}

int main(int argc, char* argv[]) {
    AnalyzeSignalForAllMasses();
    AnalyzeBackground();
    
    for (const std::pair<double, double>& p : gd.resolvedYield) {
        rbLow[p.first] = p.second + gd.boostedLowYield[p.first];
        rbHigh[p.first] = p.second + gd.boostedHighYield[p.first];
    }
    
    signResolved = CalculateSignificance(gd.resolvedYield);
    signLow = CalculateSignificance(rbLow);
    signHigh = CalculateSignificance(rbHigh);
    
    effHigh = CalculateEfficiency(gd.boostedHighCount);
    effLow = CalculateEfficiency(gd.boostedLowCount);
    effResolved = CalculateEfficiency(gd.resolvedCount);
    
    PlotYields();
    PlotSignificance();
    PlotEfficiency();
    
    AnalyzeBackgroundSlices();
    
    return 0;
}
