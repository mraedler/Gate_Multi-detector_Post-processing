#include <iostream>
#include <TFile.h>
#include "include/utils.h"
#include "include/eventSelection.h"

// Globals
bool g_verbose;
TString g_treeName;
std::vector<TString> g_fullPaths;



void processSingleRootFile(const size_t idx) {
    std::cout << "Processing: " << g_fullPaths[idx] << std::endl;

    TFile* file = getTFile(g_fullPaths[idx], "UPDATE", g_verbose);
    TTree* tree = getTTree(file, g_treeName, g_verbose);

    checkIfLeafAlreadyExists(tree, "group");

    Int_t groupID, groupEventID1, groupEventID2;
    TBranch* b0 = tree->Branch("groupID", &groupID, "groupID/I");
    TBranch* b1 = tree->Branch("groupEventID1", &groupEventID1, "groupEventID1/I");
    TBranch* b2 = tree->Branch("groupEventID2", &groupEventID2, "groupEventID2/I");

    coincidenceGrouping(tree, groupID, b0, groupEventID1, b1, groupEventID2, b2, g_verbose);

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

    // runSequentially(g_fullPaths.size(), processSingleFile);
    runInSeparateProcesses(g_fullPaths.size(), processSingleRootFile, 128);


    return 0;
}