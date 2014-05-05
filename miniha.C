#include <string>
#include <iostream>
#include <vector>

#include "DelphesNTuple.h"
#include "LocalSettings.h"
#include "HiggsHist.h"

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
  HiggsHist(inputFolder, filenames, xsections, outputFolder);
  std::cout << "Done with background" << std::endl;
}

void
AnalyzeSignalForAllMasses() {
  int masses[] = { 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000 };
  for (int i = 0; i < sizeof(masses)/sizeof(masses[0]); ++i) {
    char fileName[1024];
    sprintf(fileName, "mintree_jetsub_a-zh-triple-%dGeV.root", masses[i]);
    HiggsHist(inputFolder, fileName, outputFolder);
    std::cout << "Done with " << masses[i] << "GeV" << std::endl;
    std::cout.flush();
  }
}

int main(int argc, char* argv[]) {
  AnalyzeSignalForAllMasses();
  AnalyzeBackground();
  return 0;
}
