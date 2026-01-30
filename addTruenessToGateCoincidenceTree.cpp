#include <iostream>
#include <TFile.h>
#include "include/utils.h"
// #include "include/CASToRTools.h"
#include "include/eventSelection.h"

// Globals
bool g_verbose;
TString g_treeName;
std::vector<TString> g_fullPaths;



void processSingleRootFile(const size_t idx) {
    std::cout << "Processing: " << g_fullPaths[idx] << std::endl;

    TFile* file = getTFile(g_fullPaths[idx], "UPDATE", g_verbose);
    TTree* tree = getTTree(file, g_treeName, g_verbose);

    checkIfLeafAlreadyExists(tree, "trueness");

    // Add boolean branch for the preselection
    Bool_t trueness;
    TBranch* b = tree->Branch("trueness", &trueness, "trueness/O");

    identifyTrueEvents(tree, trueness, b);

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

    runInSeparateProcesses(g_fullPaths.size(), processSingleRootFile, 128);
    //runSequentially(g_fullPaths.size(), processSingleFile);

    return 0;
}