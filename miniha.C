#include <string>
#include <iostream>

#include "DelphesNTuple.h"
#include "LocalSettings.h"
#include "HiggsHist.h"

void
AnalyzeBackground() {
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
  return 0;
}
