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

    if (!checkIfLeafExists(tree, "scatterTest") || !checkIfLeafExists(tree, "trueness")) {
        std::cerr << "Error: coincidence attributes leafs missing. Run addCoincidenceAttributes first." << std::endl;
        std::exit(1);
    }

    // Get necessary leafs
    // todo: Can also be moved to the respective event selection functions
    TLeaf* gantryID1 = tree->GetLeaf("gantryID1");
    TLeaf* gantryID2 = tree->GetLeaf("gantryID2");
    TLeaf* rsectorID1 = tree->GetLeaf("rsectorID1");
    TLeaf* rsectorID2 = tree->GetLeaf("rsectorID2");
    TLeaf* trueness = tree->GetLeaf("trueness");
    TLeaf* scatterTest = tree->GetLeaf("scatterTest");

    // Add boolean branch for the preselection
    Bool_t preselection;
    TBranch* b = tree->Branch("preselection", &preselection, "preselection/O");

    runMinSectorDifferenceTest(tree, gantryID1, gantryID2, rsectorID1, rsectorID2, preselection, b, g_minSectorDifference, g_verbose);
    runScatterTest(tree, scatterTest, gantryID1, gantryID2, trueness, preselection, b, g_verbose);

    // Set up new file containing only the preselected data
    TString newFullPath = g_fullPaths[idx];
    newFullPath.Remove(newFullPath.Last('.'));
    newFullPath += "_new.root";
    TFile* newFile = TFile::Open(newFullPath, "RECREATE");
    //newFile->cd();  // Ensure the new file is linked
    TTree* newTree = tree->CloneTree(0);
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t ii = 0; ii < nEntries; ++ii) {
        tree->GetEntry(ii);
        if (!preselection) continue;
        newTree->Fill();
    }

    std::cout << "Writing: " << newFullPath << std::endl;

    newTree->Write();
    newFile->Close();
    file->Close();

    std::remove(g_fullPaths[idx]);
    std::rename(newFullPath, g_fullPaths[idx]);  // rename new file
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
    // todo: If you change the LUT here, also change in the runMinSectorDifferenceTest (fix this)
    // g_lut = readLutBinary("/data/local1/raedler/J-PET/CASToR/castor/config/scanner/TB_J-PET_7th_gen_brain_insert_dz_1_mm.lut");
    // g_lut = readLutBinary("/data/local1/raedler/J-PET/CASToR/castor/config/scanner/TB_J-PET_7th_gen_brain_insert_WHR_4_18_1_mm.lut");
    g_lut = readLutBinary("/data/local1/raedler/J-PET/CASToR/castor/config/scanner/TB_J-PET_7th_gen_brain_insert_WHR_6_30_1_mm.lut");

    // runSequentially(g_fullPaths.size(), processSingleFile);
    runInSeparateProcesses(g_fullPaths.size(), processSingleFile, 128);

    return 0;
}