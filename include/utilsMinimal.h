#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <TString.h>
#include <TTree.h>
#include <TLeaf.h>

class TFile;

void listAvailableTrees(const TFile* file);
void listAvailableBranches(TTree* tree);
void listAvailableLeaves(TTree* tree);
bool checkIfLeafExists(TTree* tree, const std::string& checkLeafName);

std::vector<TString> getListOfRootFilePaths(const char* dirPath, bool verbose = false);

using JobFn = void(*)(size_t);
void runInSeparateProcesses(size_t nJobs, JobFn job, int maxProcesses = 128);
void runSequentially(size_t nJobs, JobFn job);

TFile* getTFile(const TString& fullPath, const char* mode, bool verbose = false);
TTree* getTTree(TFile* file, const TString& treeName, bool verbose = false);

#endif //UTILS_H
