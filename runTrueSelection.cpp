#include <iomanip>
#include <iostream>
#include <TFile.h>
#include "include/utilsMinimal.h"
#include "include/eventSelection.h"

// Globals
bool g_verbose;
TString g_treeName;
std::vector<TString> g_fullPaths;



bool checkLeaf(TTree* tree, const char* leafName) {
    bool exists = true;
    if (!checkIfLeafExists(tree, leafName)) {
        std::cout << "Missing leaf '" << leafName << "'." << std::endl;
        exists = false;
    }
    return exists;
}



void processSingleFile(const size_t idx) {
    std::cout << "Processing: " << g_fullPaths[idx] << std::endl;

    TFile* file = getTFile(g_fullPaths[idx], "UPDATE", g_verbose);
    TTree* tree = getTTree(file, g_treeName, g_verbose);

    bool allNecessaryLeafsExist = (
        checkLeaf(tree, "eventID1") &&
        checkLeaf(tree, "eventID2") &&
        checkLeaf(tree, "comptonCrystal1") &&
        checkLeaf(tree, "comptonCrystal2") &&
        checkLeaf(tree, "RayleighCrystal1") &&
        checkLeaf(tree, "RayleighCrystal2") &&
        checkLeaf(tree, "comptonPhantom1") &&
        checkLeaf(tree, "comptonPhantom2") &&
        checkLeaf(tree, "RayleighPhantom1") &&
        checkLeaf(tree, "RayleighPhantom2") &&
        checkLeaf(tree, "sourcePosX1") &&
        checkLeaf(tree, "sourcePosY1") &&
        checkLeaf(tree, "sourcePosZ1"));

    if (!allNecessaryLeafsExist) {
        std::cout << "Missing leaf(s). Skipping '" << g_fullPaths[idx] << "'." << std::endl;
        return;
    }

    // Add boolean branch indicating the "trueness"
    Bool_t trueness;
    TBranch* b = tree->Branch("trueness", &trueness, "trueness/O");

    // Fill the branch
    identifyTrueEvents(tree, trueness, b);

    // Deactivate all branches except for the source position and the added trueness branch
    tree->SetBranchStatus("*", false);
    tree->SetBranchStatus("sourcePosX1", true);
    tree->SetBranchStatus("sourcePosY1", true);
    tree->SetBranchStatus("sourcePosZ1", true);
    // tree->SetBranchStatus("sourcePosX2", true);
    // tree->SetBranchStatus("sourcePosY2", true);
    // tree->SetBranchStatus("sourcePosZ2", true);

    // Set up new file containing only the preselected data
    TString newFullPath = g_fullPaths[idx];
    newFullPath.Remove(newFullPath.Last('.'));
    newFullPath += "_True.root";
    TFile* newFile = TFile::Open(newFullPath, "RECREATE");
    //newFile->cd();  // Ensure the new file is linked
    TTree* newTree = tree->CloneTree(0);

    //
    tree->SetBranchStatus("trueness", true);

    Long64_t nEntries = tree->GetEntries();
    Long64_t progressStep = nEntries / 10;

    for (Long64_t ii = 0; ii < nEntries; ++ii) {
        tree->GetEntry(ii);

        // Print progress (adding nEntries / 2 corresponds to rounding in integer division)
        if (ii % progressStep == 0) {std::cout << (ii * 100 + nEntries / 2) / nEntries << " %" << std::endl;}

        if (!trueness) {continue;}

        newTree->Fill();
    }


    std::cout << "Writing: " << newFullPath << std::endl;

    newTree->Write("", TObject::kOverwrite);
    newFile->Close();
    file->Close();

    std::remove(g_fullPaths[idx]);
    // std::rename(newFullPath, g_fullPaths[idx]);  // rename new file
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

    // runSequentially(g_fullPaths.size(), processSingleFile);
    runInSeparateProcesses(g_fullPaths.size(), processSingleFile, 128);

    return 0;
}