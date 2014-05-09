#ifndef _HIGGS_HIST_H_
#define _HIGGS_HIST_H_

void HiggsHist(const char* inputFolder, const char* inputFile, const char* outputFolder, double xsection, double mass, double eventCount);
void HiggsHist(const char* inputFolder, std::vector<const char*>& inputFiles, std::vector<double>& crossSections, const char* outputFolder, double mass, double eventCountPerFile);

#endif
