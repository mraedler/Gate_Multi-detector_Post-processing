#include <iostream>
#include <TFile.h>
#include "include/utils.h"
#include "include/CASToRTools.h"

// Globals
bool g_verbose;
TString g_treeName;
std::vector<TString> g_fullPaths;
std::vector<LutEntry> g_lut;



void processSingleRootFile(const size_t idx) {
    // todo: avoid writing in the castorID vectors; write directly into the new branches
    std::cout << "Processing: " << g_fullPaths[idx] << std::endl;

    TFile* file = getTFile(g_fullPaths[idx], "UPDATE", g_verbose);
    TTree* tree = getTTree(file, g_treeName, g_verbose);

    checkIfLeafAlreadyExists(tree, "castorID");

    // For the indexing
    TLeaf* gantryID1 = tree->GetLeaf("gantryID1");
    TLeaf* rsectorID1 = tree->GetLeaf("rsectorID1");
    TLeaf* crystalID1 = tree->GetLeaf("crystalID1");
    TLeaf* layerID1 = tree->GetLeaf("layerID1");

    // For the consistency check
    // TLeaf* globalPosX1 = tree->GetLeaf("globalPosX1");
    // TLeaf* globalPosY1 = tree->GetLeaf("globalPosY1");
    // TLeaf* globalPosZ1 = tree->GetLeaf("globalPosZ1");

    std::vector<int> vectorCastorID1 = getCastorID(tree, gantryID1, rsectorID1, crystalID1, layerID1);
    //checkCastorID(tree, gantryID1, globalPosX1, globalPosY1, globalPosZ1, vectorCastorID1, g_lut);

    TLeaf* gantryID2 = tree->GetLeaf("gantryID2");
    TLeaf* rsectorID2 = tree->GetLeaf("rsectorID2");
    TLeaf* crystalID2 = tree->GetLeaf("crystalID2");
    TLeaf* layerID2 = tree->GetLeaf("layerID2");

    // TLeaf* globalPosX2 = tree->GetLeaf("globalPosX2");
    // TLeaf* globalPosY2 = tree->GetLeaf("globalPosY2");
    // TLeaf* globalPosZ2 = tree->GetLeaf("globalPosZ2");

    std::vector<int> vectorCastorID2 = getCastorID(tree, gantryID2, rsectorID2, crystalID2, layerID2);
    //checkCastorID(tree, gantryID2, globalPosX2, globalPosY2, globalPosZ2, vectorCastorID2, g_lut);

    // clone structure + data
    //TTree* newTree = tree->CloneTree(0);

    Int_t castorID1, castorID2;
    TBranch* b1 = tree->Branch("castorID1", &castorID1, "castorID1/I");
    TBranch* b2 = tree->Branch("castorID2", &castorID2, "castorID2/I");

    Long64_t nEntries = tree->GetEntries();
    if ((vectorCastorID1.size() != nEntries) || (vectorCastorID2.size() != nEntries)) {
        std::cout << "Error: vector size mismatch." << std::endl;
        file->Close();
        std::exit(1);
    }

    for (Long64_t ii = 0; ii < nEntries; ++ii) {
        castorID1 = vectorCastorID1[ii];
        castorID2 = vectorCastorID2[ii];
        b1->Fill();
        b2->Fill();

        // Progress bar
        // if (ii % (nEntries / 100) == 0) {  // update every 1%
            // int progress = (ii * 100) / nEntries;
            // std::cout << "\rProgress: " << progress << "%";
            // std::cout.flush();
        // }
    }

    //sleep(100);

    std::cout << "Writing: " << g_fullPaths[idx] << std::endl;

    file->cd();
    // tree->Write(g_treeName, TObject::kOverwrite);  // Overwrite
    // tree->Write("", TObject::kOverwrite);  // Overwrite
    tree->Write("", TObject::kWriteDelete);
    file->Close();
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
    g_lut = readLutBinary("/data/local1/raedler/J-PET/CASToR/castor/config/scanner/TB_J-PET_7th_gen_brain_insert_dz_1_mm.lut");

    //runSequentially(g_fullPaths.size(), processSingleFile);
    runInSeparateProcesses(g_fullPaths.size(), processSingleRootFile, 128);

    return 0;
}