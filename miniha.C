#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "DelphesNTuple.h"
#include "LocalSettings.h"
#include "HiggsHist.h"
#include "HelperClasses.h"

void
AnalyzeBackground() {
  std::vector<const char*> filenames;
  std::vector<double> xsections;
  for (int i = 1; i <= 8; ++i) {
    char buffer[2048];
    sprintf(buffer, "mintree_jetsub_back-triple-%dGeV.root", i);
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
  HiggsHist(inputFolder, filenames, xsections, outputFolder, 1, 50000);
  std::cout << "Done with background" << std::endl;
}

void
AnalyzeSignalForAllMasses() {
  int masses[] = { 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000 };
  double xsections[] = { 0.405, 0.493, 0.0473, 0.0114, 0.00431, 0.00205, 0.00106, 0.000577, 0.000326, 0.000189 };
  for (int i = 0; i < sizeof(masses)/sizeof(masses[0]); ++i) {
    char fileName[1024];
    sprintf(fileName, "mintree_jetsub_a-zh-triple-%dGeV.root", masses[i]);
    HiggsHist(inputFolder, fileName, outputFolder, xsections[i], masses[i], 50000);
    std::cout << "Done with " << masses[i] << "GeV" << std::endl;
    std::cout.flush();
  }
}

void
PlotYields(const std::string& name, std::map<double, double> yields) {
  double miny = std::min_element(yields.begin(), yields.end(), [](std::pair<double, double> p, std::pair<double, double> q) { return p.second < q.second; })->second;
  double maxy = std::max_element(yields.begin(), yields.end(), [](std::pair<double, double> p, std::pair<double, double> q) { return p.second < q.second; })->second;

  TH2F* hist = new TH2F(name.c_str(), name.c_str(), 100, 0, 1000, 100, miny * 0.9, maxy * 1.1);
  for (const std::pair<double, double>& p : yields) {
    hist->Fill(p.first, p.second);
  }
}

int main(int argc, char* argv[]) {
  AnalyzeSignalForAllMasses();
  AnalyzeBackground();

  TFile* yieldsFile = new TFile(outputFolder + "yeilds.root", "recreate");
  PlotYields("resolved", gd.resolvedYield);
  PlotYields("boosted", gd.boostedYield);
  yieldsFile->Write();
  yieldsFile->Close();

  return 0;
}
