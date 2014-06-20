#include <TString.h>

TString outputFolder = "/Research/outputs/";
TString inputFolder = "/Research/outputs/";
TString libDelphesPath = "/phys/groups/tev/scratch3/users/milenl/MadGraph5/Delphes/libDelphes";

//int amasses[] = { 3000 };
int amasses[] = { 300, 350, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3000 };
double xsections[] = { 0.493, 0.0473, 0.0114, 0.00431, 0.00205, 0.00106, 0.000577, 0.000326, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189, 0.000189};

int amassesCount = sizeof(amasses)/sizeof(amasses[0]);
