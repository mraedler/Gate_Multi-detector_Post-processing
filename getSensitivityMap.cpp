#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TH3D.h>
#include <TSystem.h>
#include "TROOT.h"
#include "include/utils.h"
// #include "include/eventSelection.h"

// Globals
bool g_verbose;
TString g_treeName;
std::vector<TString> g_fullPaths;
TString g_outputPath;

std::vector<int> g_nVoxels;
std::vector<double> g_mapCenter;
std::vector<double> g_mapHalfSize;


void processSingleFile(const size_t idx) {
    std::cout << "Processing: " << g_fullPaths[idx] << std::endl;

    TFile* file = getTFile(g_fullPaths[idx], "UPDATE", g_verbose);
    TTree* tree = getTTree(file, g_treeName, g_verbose);
    // listAvailableBranches(tree);

    TString fileName = gSystem->BaseName(g_fullPaths[idx]);
    fileName.Remove(fileName.Last('.'));

    // Allocate the output histograms
    TH3D* hist = new TH3D(fileName + "_hist", "Counts",
        g_nVoxels[0], g_mapCenter[0] - g_mapHalfSize[0], g_mapCenter[0] + g_mapHalfSize[0],
        g_nVoxels[1], g_mapCenter[1] - g_mapHalfSize[1], g_mapCenter[1] + g_mapHalfSize[1],
        g_nVoxels[2], g_mapCenter[2] - g_mapHalfSize[2], g_mapCenter[2] + g_mapHalfSize[2]);

    //
    TLeaf* energy1 = tree->GetLeaf("energy1");
    TLeaf* energy2 = tree->GetLeaf("energy2");
    TLeaf* eventID1 = tree->GetLeaf("eventID1");
    TLeaf* eventID2 = tree->GetLeaf("eventID2");
    TLeaf* comptonCrystal1 = tree->GetLeaf("comptonCrystal1");
    TLeaf* comptonCrystal2 = tree->GetLeaf("comptonCrystal2");
    TLeaf* rayleighCrystal1 = tree->GetLeaf("RayleighCrystal1");
    TLeaf* rayleighCrystal2 = tree->GetLeaf("RayleighCrystal2");

    TLeaf* sourcePosX1 = tree->GetLeaf("sourcePosX1");
    TLeaf* sourcePosY1 = tree->GetLeaf("sourcePosY1");
    TLeaf* sourcePosZ1 = tree->GetLeaf("sourcePosZ1");
    // TLeaf* sourcePosX2 = tree->GetLeaf("sourcePosX2");
    // TLeaf* sourcePosY2 = tree->GetLeaf("sourcePosY2");
    // TLeaf* sourcePosZ2 = tree->GetLeaf("sourcePosZ2");

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t ii = 0; ii < nEntries; ++ii) {
        tree->GetEntry(ii);

        bool aboveEnergyThreshold = (energy1->GetValue() > .2) && (energy2->GetValue() > .2);  // both above 200 keV

        bool sameEvent = eventID1->GetValue() == eventID2->GetValue();
        bool firstCompton = (comptonCrystal1->GetValue() == 1) && (comptonCrystal2->GetValue() == 1);
        bool zeroRayleigh = (rayleighCrystal1->GetValue() == 0) && (rayleighCrystal2->GetValue() == 0);
        bool trueCoincidence = sameEvent && firstCompton && zeroRayleigh;

        if (trueCoincidence && aboveEnergyThreshold) {
            hist->Fill(sourcePosX1->GetValue(), sourcePosY1->GetValue(), sourcePosZ1->GetValue());
        }
    }

    hist->SaveAs(g_outputPath + fileName + "_hist.root");

    delete hist;

    file->Close();
}



TH3D* getTH3D(const TString& fileName) {
    TFile* file = TFile::Open(g_outputPath +  fileName + ".root", "READ");
    TH3D* hist = (TH3D*)file->Get(fileName);
    hist->SetDirectory(nullptr);  // detach from file
    file->Close();
    delete file;
    gSystem->Unlink(g_outputPath +  fileName + ".root");
    return hist;
}



void mergeHists() {
    TString firstFileName = gSystem->BaseName(g_fullPaths[0]);
    firstFileName.Remove(firstFileName.Last('.'));

    TH3D* firstHist = getTH3D(firstFileName + "_hist");

    for (unsigned int ii = 1; ii < g_fullPaths.size(); ii++) {
        TString fileName = gSystem->BaseName(g_fullPaths[ii]);
        fileName.Remove(fileName.Last('.'));

        TH3D* hist = getTH3D(fileName + "_hist");
        firstHist->Add(hist);
        delete hist;
    }
    firstHist->SetName("TH3");
    firstHist->SaveAs(g_outputPath + "merged_hist.root");
    delete firstHist;
}



int main(int argc, char* argv[]) {
    // Check for the required input
    if (argc < 2) {
        printf("Error: no path provided.\n");
        return 1;
    }
    const char *path = argv[1];
    //printf("Path: %s\n", path);

    // Set globals
    g_verbose = false;
    g_treeName = "MergedCoincidences";
    g_fullPaths = getListOfRootFilePaths(path, g_verbose);

    // Adapt the same folder structure (last two folders) from the input
    g_outputPath = "/net/people/plgrid/plgraedler/SensitivityMap/Output/";  // needs to have trailing "/"

    g_nVoxels = {800, 800, 2540};  // 1×1×1 mm spacing
    g_mapCenter = {0., 0., 0.};  // mm
    g_mapHalfSize = {400., 400., 1270.};  // mm

    runSequentially(g_fullPaths.size(), processSingleFile);
    // runInSeparateProcesses(g_fullPaths.size(), processSingleFile, 128);

    std::exit(1);

    mergeHists();

    return 0;
}