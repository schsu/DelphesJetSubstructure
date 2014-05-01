#include <vector>
#include <thread>
#include <cstdlib>
#include <unistd.h>

#include "LocalSettings.h"
#include "HiggsAnalysis.h"
#include "JetSubstructure.h"

using namespace std;

void exec_cmd_no_wait(const std::string& cmd, std::vector<std::string>& params) {
  pid_t procId = fork();
  // if prodId == 0 -> we forked, this is the new process
  if (procId == 0) {
    char** args = new char*[params.size() + 2];
    args[0] = strdup(cmd.c_str());
    args[params.size() + 1] = NULL;
    for (size_t i = 0; i < params.size(); ++i) {
      args[i+1] = strdup(params[i].c_str());
      execvp(cmd.c_str(), args);
    }
  }
}

void RunFullBackground() {
  TString inputFile = "background-full-21-cone06.root";
  HiggsAnalysis(inputFolder, inputFile, outputFolder, 0.083, 0.6);
  inputFile = "background-full-21-cone10.root";
  HiggsAnalysis(inputFolder, inputFile, outputFolder, 0.083, 1.0);
}

void RunHeavyHiggsThread(TString inputFile) {
  JetSubstructure(inputFolder, inputFile, outputFolder, 1.0, 0.0);
}

void RunHeavyHiggs(TString coneSize, float fcs) {
  TString baseName = "a-zh-triple-";
  TString masses[] = {
    "250GeV",
    "300GeV",
    "350GeV",
    "400GeV",
    "500GeV",
    "600GeV",
    "700GeV",
    "800GeV",
    "900GeV",
    "1000GeV",
    // "1100GeV",
    // "1200GeV",
    // "1300GeV",
    // "1400GeV",
    // "1500GeV",
  };

  for (int i = 0; i < sizeof(masses)/sizeof(masses[0]); ++i) {
    TString inputFile = baseName + masses[i] + ".root";
    std::vector<std::string> params;
    params.push_back(std::string((const char*)inputFile));
    exec_cmd_no_wait("./ha", params);
  }
}

void RunBackground() {
  TString baseName = "back-triple-";
  TString masses[] = {
    "1GeV",
    "2GeV",
    "3GeV",
    "4GeV",
    "5GeV",
    "6GeV",
    "7GeV",
    "8GeV",
  };

  for (int i = 0; i < sizeof(masses)/sizeof(masses[0]); ++i) {
    TString inputFile = baseName + masses[i] + ".root";
    std::vector<std::string> params;
    params.push_back(std::string((const char*)inputFile));
    exec_cmd_no_wait("./ha", params);
  }
}

void RunBackground10() {
  vector<TString> inputs;
  inputs.push_back("background-1-cone10.root");
  inputs.push_back("background-2-cone10.root");
  inputs.push_back("background-3-cone10.root");
  inputs.push_back("background-4-cone10.root");
  inputs.push_back("background-5-cone10.root");
  inputs.push_back("background-6-cone10.root");
  inputs.push_back("background-7-cone10.root");
  inputs.push_back("background-8-cone10.root");
  
  vector<float> xsections;
  xsections.push_back(14.59);
  xsections.push_back(0.7155);
  xsections.push_back(0.0668);
  xsections.push_back(0.01246);
  xsections.push_back(0.003444);
  xsections.push_back(0.001215);
  xsections.push_back(0.0007277);
  xsections.push_back(0.0001713);

  BachgroundAnalysis(inputFolder, inputs, xsections, outputFolder, 1.0);
}

void RunTruthMass(TString mass, float m) {
  TString inputFile = "a-zh-h125-" + mass + "GeV-cone06.root";
  HiggsTruthAnalysis(inputFolder, inputFile, outputFolder, m);
}

void RunTruth() {
  RunTruthMass("300", 300.0);
  RunTruthMass("350", 350.0);
  RunTruthMass("500", 500.0);
  RunTruthMass("800", 800.0);
  RunTruthMass("1000", 1000.0);
}

void Run() {
  //gSystem->Load(libDelphesPath);
  //gSystem->CompileMacro("HelperClasses.cxx");
  //gSystem->CompileMacro("HiggsAnalysis.C");
  //gSystem->CompileMacro("JetSubstructure.C");
  //gSystem->CompileMacro("Qjets.C");
  //RunBackground06();
  //RunHeavyHiggsExtraCuts10();
  //RunHeavyHiggsExtraCuts06();
  //RunHeavyHiggs("", 1.0);
  RunBackground();
  //RunHeavyHiggs10();
  //RunTruth();
  //RunFullBackground();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    Run();
  } else {
    cout << argv[1] << endl;
    RunHeavyHiggsThread(argv[1]);
  }

  return 0;
}
