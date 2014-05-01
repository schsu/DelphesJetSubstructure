#ifndef _HIGGS_ANALYSIS_H_
#define _HIGGS_ANALYSIS_H_

#include <vector>
using namespace std;

int HiggsAnalysisExtraCuts(const char* inputFolder, const char* inputFile, const char* outputFolder, float xsection, float amass, float coneSize);
int HiggsAnalysis(const char* inputFolder, const char* inputFile, const char* outputFolder, float xsection, float coneSize);
int HiggsTruthAnalysis(const char* inputFolder, const char* inputFile, const char* outputFolder, float heavyHiggsMass);
int HiggsAnalysisNoCuts(const char* inputFolder, const char* intputFile);
int HiggsPlusBackground(const char* inputFolder, const char* inputFile, const char* inputBackgroundFolder, const char* inputBackgroundFile);
int BachgroundAnalysis(const char* inputFolder, vector<TString>& inputFile, vector<float>& xsections, const char* outputFolder, float coneSize);

#endif
