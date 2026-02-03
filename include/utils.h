#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <TString.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1D.h>

class TFile;
// TTree;
//class TString;

void listAvailableTrees(const TFile* file);
void listAvailableBranches(TTree* tree);
void listAvailableLeaves(TTree* tree);
void checkIfLeafAlreadyExists(TTree* tree, const std::string& checkLeafName);

std::vector<TString> getListOfRootFilePaths(const char* dirPath, bool verbose = false);

using JobFn = void(*)(size_t);
void runInSeparateProcesses(size_t nJobs, JobFn job, int maxProcesses = 128);
void runSequentially(size_t nJobs, JobFn job);

TFile* getTFile(const TString& fullPath, const char* mode, bool verbose = false);
TTree* getTTree(TFile* file, const TString& treeName, bool verbose = false);

template <typename T>
TH1D* getHistogram(std::vector<T>& arr, const float lowerEdge, const float binWidth, const int nBins, const char* name, const char* labeling) {

    TH1D* hist = new TH1D(name, labeling, nBins, lowerEdge, lowerEdge + nBins * binWidth);
    for (const T& a : arr)
        hist->Fill(a);

    return hist;
}

std::vector<TString> splitPath(const TString& path);

#endif //UTILS_H
