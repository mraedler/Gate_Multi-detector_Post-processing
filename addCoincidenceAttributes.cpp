#include <iostream>
#include <TFile.h>
#include "include/utils.h"
#include "include/eventSelection.h"
#include "include/CASToRTools.h"

// Globals
bool g_verbose;
TString g_treeName;
std::vector<TString> g_fullPaths;
TString g_lutName;
std::vector<LutEntry> g_lut;



void processSingleFile(const size_t idx) {
    std::cout << "Processing: " << g_fullPaths[idx] << std::endl;

    TFile* file = getTFile(g_fullPaths[idx], "UPDATE", g_verbose);
    TTree* tree = getTTree(file, g_treeName, g_verbose);

    // Flag if the event is true
    if (!checkIfLeafExists(tree, "trueness")) {
        // Add boolean branch for the preselection
        Bool_t trueness;
        TBranch* b = tree->Branch("trueness", &trueness, "trueness/O");

        // Fill the branch
        identifyTrueEvents(tree, trueness, b);
    }

    // Index for CASToR LUT
    Int_t castorID1, castorID2;
    if (!checkIfLeafExists(tree, "castorID")) {
        TBranch* b1 = tree->Branch("castorID1", &castorID1, "castorID1/I");
        TBranch* b2 = tree->Branch("castorID2", &castorID2, "castorID2/I");

        // For the indexing
        TLeaf* gantryID1 = tree->GetLeaf("gantryID1");
        TLeaf* rsectorID1 = tree->GetLeaf("rsectorID1");
        TLeaf* crystalID1 = tree->GetLeaf("crystalID1");
        TLeaf* layerID1 = tree->GetLeaf("layerID1");

        // // For the consistency check
        // TLeaf* globalPosX1 = tree->GetLeaf("globalPosX1");
        // TLeaf* globalPosY1 = tree->GetLeaf("globalPosY1");
        // TLeaf* globalPosZ1 = tree->GetLeaf("globalPosZ1");

        setCastorID(tree, gantryID1, rsectorID1, crystalID1, layerID1, castorID1, b1, g_lutName);
        // checkCastorID(tree, gantryID1, globalPosX1, globalPosY1, globalPosZ1, castorID1, g_lut);

        TLeaf* gantryID2 = tree->GetLeaf("gantryID2");
        TLeaf* rsectorID2 = tree->GetLeaf("rsectorID2");
        TLeaf* crystalID2 = tree->GetLeaf("crystalID2");
        TLeaf* layerID2 = tree->GetLeaf("layerID2");

        // TLeaf* globalPosX2 = tree->GetLeaf("globalPosX2");
        // TLeaf* globalPosY2 = tree->GetLeaf("globalPosY2");
        // TLeaf* globalPosZ2 = tree->GetLeaf("globalPosZ2");

        setCastorID(tree, gantryID2, rsectorID2, crystalID2, layerID2, castorID2, b2, g_lutName);
        // checkCastorID(tree, gantryID2, globalPosX2, globalPosY2, globalPosZ2, castorID2, g_lut);
    }
    else {
        tree->SetBranchAddress("castorID1", &castorID1);
        tree->SetBranchAddress("castorID2", &castorID2);
    }

    // Scatter test
    if (!checkIfLeafExists(tree, "scatterTest")) {
        Double_t scatterTest;
        TBranch* b = tree->Branch("scatterTest", &scatterTest, "scatterTest/D");

        setScatterTest(tree, castorID1, castorID2, g_lut, scatterTest, b, g_verbose);
    }

    std::cout << "Writing: " << g_fullPaths[idx] << std::endl;

    file->cd();
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

    g_lutName =  "TB_J-PET_7th_gen_brain_insert_WHR_4_18_1_mm";
    std::set<TString> lutOptions = {"TB_J-PET_7th_gen_brain_insert_WHR_4_18_1_mm", "TB_J-PET_7th_gen_brain_insert_WHR_6_30_1_mm"};

    bool allArgsProvidedAndUsed = argc >= 3;  // todo: increment if more optional args are added

    // Parse optional arguments
    for (int ii = 2; ii < argc; ++ii) {
        std::string arg = argv[ii];

        bool argRecognized = false;
        setArgument("lut", arg, lutOptions, g_lutName, argRecognized, allArgsProvidedAndUsed);

        if (!argRecognized) {std::cout << "Warning: unknown argument '" << arg << "' ignored.\n";}
    }

    if (!allArgsProvidedAndUsed) {
        std::cout << "\nWarning: not all optional arguments were provided or valid.\n";
        std::cout << "Using the following (default) parameters:\n";
        std::cout << "  lut = " << g_lutName << "\n\n";

        std::cout << "Press ENTER to continue or Ctrl+C to abort...";
        std::cin.get();  // wait for Enter
    }

    g_lut = readLutBinary("/data/local1/raedler/J-PET/CASToR/castor/config/scanner/" + g_lutName + ".lut");

    // runSequentially(g_fullPaths.size(), processSingleFile);
    runInSeparateProcesses(g_fullPaths.size(), processSingleFile, 128);

    return 0;
}