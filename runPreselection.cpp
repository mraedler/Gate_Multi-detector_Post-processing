#include <iostream>
#include <TFile.h>
#include "include/utils.h"
#include "include/CASToRTools.h"
#include "include/eventSelection.h"

// Globals
bool g_verbose;
int g_minSectorDifference;
TString g_treeName;
std::vector<TString> g_fullPaths;
std::vector<LutEntry> g_lut;



void processSingleFile(const size_t idx) {
    std::cout << "Processing: " << g_fullPaths[idx] << std::endl;

    TFile* file = getTFile(g_fullPaths[idx], "UPDATE", g_verbose);
    TTree* tree = getTTree(file, g_treeName, g_verbose);

    //
    TString newFullPath = g_fullPaths[idx];
    newFullPath.Remove(newFullPath.Last('.'));
    newFullPath += "_new.root";
    TFile* newFile = TFile::Open(newFullPath, "RECREATE");
    //newFile->cd();  // Ensure the new file is linked
    TTree* newTree = tree->CloneTree(0);

    //
    TLeaf* gantryID1 = tree->GetLeaf("gantryID1");
    TLeaf* gantryID2 = tree->GetLeaf("gantryID2");
    TLeaf* rsectorID1 = tree->GetLeaf("rsectorID1");
    TLeaf* rsectorID2 = tree->GetLeaf("rsectorID2");
    TLeaf* castorID1 = tree->GetLeaf("castorID1");
    TLeaf* castorID2 = tree->GetLeaf("castorID2");
    TLeaf* time1 = tree->GetLeaf("time1");
    TLeaf* time2 = tree->GetLeaf("time2");
    TLeaf* trueness = tree->GetLeaf("trueness");

    // Add boolean branch for the preselection
    Bool_t preselection;
    TBranch* b = tree->Branch("preselection", &preselection, "preselection/O");

    runMinSectorDifferenceTest(tree, gantryID1, gantryID2, rsectorID1, rsectorID2, preselection, b, g_minSectorDifference, false);
    runScatterTest(tree, time1, time2, castorID1, castorID2, gantryID1, gantryID2, trueness, g_lut, preselection, b, false);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t ii = 0; ii < nEntries; ++ii) {
        tree->GetEntry(ii);
        if (!preselection) continue;
        newTree->Fill();
    }

    newTree->Write();
    newFile->Close();
    file->Close();

    // #include <cstdio>
    std::remove(g_fullPaths[idx]);
    std::rename(newFullPath, g_fullPaths[idx]);  // rename new file

    // std::exit(1);
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
    g_minSectorDifference = 2;
    g_treeName = "MergedCoincidences";
    g_fullPaths = getListOfRootFilePaths(path, g_verbose);
    g_lut = readLutBinary("/data/local1/raedler/J-PET/CASToR/castor/config/scanner/TB_J-PET_7th_gen_brain_insert_dz_1_mm.lut");

    runInSeparateProcesses(g_fullPaths.size(), processSingleFile, 128);
    //runSequentially(g_fullPaths.size(), processSingleFile);

    return 0;
}