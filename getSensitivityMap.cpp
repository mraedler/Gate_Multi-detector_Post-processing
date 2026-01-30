#include <iostream>
#include <TFile.h>
#include <TH3D.h>
#include <TSystem.h>
#include "include/utils.h"
#include "include/eventSelection.h"

// Globals
bool g_verbose;
TString g_treeName;
std::vector<TString> g_fullPaths;
TString g_outputPath;
std::vector<double> g_mapCenter;
std::vector<double> g_mapHalfSize;
std::vector<int> g_nVoxels;



void processSingleFile(const size_t idx) {
    std::cout << "Processing: " << g_fullPaths[idx] << std::endl;

    TFile* file = getTFile(g_fullPaths[idx], "READ", g_verbose);
    TTree* tree = getTTree(file, g_treeName, g_verbose);

    // Filename for output
    std::vector<TString> fullPathSeparated = separatePath(g_fullPaths[idx]);
    size_t nSep = fullPathSeparated.size();
    TString baseName = fullPathSeparated[nSep - 1].Remove(fullPathSeparated[nSep - 1].Last('.'));
    TString outputPath = g_outputPath + fullPathSeparated[nSep - 3] + "/" + fullPathSeparated[nSep - 2] + "/";
    std::cout << outputPath << std::endl;
    if (gSystem->AccessPathName(outputPath)) {gSystem->mkdir(outputPath, true);}

    // Allocate the output arrays
    TH3D* h_TBTB = new TH3D(baseName + "_TBTB", "Counts",
                       g_nVoxels[0], g_mapCenter[0] - g_mapHalfSize[0], g_mapCenter[0] + g_mapHalfSize[0],
                       g_nVoxels[1], g_mapCenter[1] - g_mapHalfSize[1], g_mapCenter[1] + g_mapHalfSize[1],
                       g_nVoxels[2], g_mapCenter[2] - g_mapHalfSize[2], g_mapCenter[2] + g_mapHalfSize[2]);

    TH3D* h_TBBI = new TH3D(baseName + "_TBBI", "Counts",
                       g_nVoxels[0], g_mapCenter[0] - g_mapHalfSize[0], g_mapCenter[0] + g_mapHalfSize[0],
                       g_nVoxels[1], g_mapCenter[1] - g_mapHalfSize[1], g_mapCenter[1] + g_mapHalfSize[1],
                       g_nVoxels[2], g_mapCenter[2] - g_mapHalfSize[2], g_mapCenter[2] + g_mapHalfSize[2]);

    TH3D* h_BIBI = new TH3D(baseName + "_BIBI", "Counts",
                       g_nVoxels[0], g_mapCenter[0] - g_mapHalfSize[0], g_mapCenter[0] + g_mapHalfSize[0],
                       g_nVoxels[1], g_mapCenter[1] - g_mapHalfSize[1], g_mapCenter[1] + g_mapHalfSize[1],
                       g_nVoxels[2], g_mapCenter[2] - g_mapHalfSize[2], g_mapCenter[2] + g_mapHalfSize[2]);

    TLeaf* trueness = tree->GetLeaf("trueness");

    TLeaf* energy1 = tree->GetLeaf("energy1");
    TLeaf* energy2 = tree->GetLeaf("energy2");

    TLeaf* gantryID1 = tree->GetLeaf("gantryID1");
    TLeaf* gantryID2 = tree->GetLeaf("gantryID2");

    TLeaf* sourcePosX1 = tree->GetLeaf("sourcePosX1");
    TLeaf* sourcePosY1 = tree->GetLeaf("sourcePosY1");
    TLeaf* sourcePosZ1 = tree->GetLeaf("sourcePosZ1");
    //TLeaf* sourcePosX2 = tree->GetLeaf("sourcePosX2");
    //TLeaf* sourcePosY2 = tree->GetLeaf("sourcePosY2");
    //TLeaf* sourcePosZ2 = tree->GetLeaf("sourcePosZ2");

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t ii = 0; ii < nEntries; ++ii) {
        tree->GetEntry(ii);
        bool aboveEnergyThreshold = (energy1->GetValue() > .2) && (energy2->GetValue() > .2);  // both above 200 keV

        if (trueness->GetValue() && aboveEnergyThreshold) {
            if ((gantryID1->GetValue() < 2) && (gantryID2->GetValue() < 2)) {
                h_TBTB->Fill(sourcePosX1->GetValue(), sourcePosY1->GetValue(), sourcePosZ1->GetValue());
            } else if (((gantryID1->GetValue() < 2) && (gantryID2->GetValue() == 2)) || ((gantryID1->GetValue() == 2) && (gantryID2->GetValue() < 2))) {
                h_TBBI->Fill(sourcePosX1->GetValue(), sourcePosY1->GetValue(), sourcePosZ1->GetValue());
            } else if ((gantryID1->GetValue() == 2) && (gantryID2->GetValue() == 2)) {
                h_BIBI->Fill(sourcePosX1->GetValue(), sourcePosY1->GetValue(), sourcePosZ1->GetValue());
            }
        }
    }

    h_TBTB->SaveAs(outputPath + baseName + "_TBTB.root");
    h_TBBI->SaveAs(outputPath + baseName + "_TBBI.root");
    h_BIBI->SaveAs(outputPath + baseName + "_BIBI.root");

    delete h_TBTB;

    file->Close();
    std::exit(1);
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
    g_verbose = true;
    g_treeName = "MergedCoincidences";
    g_fullPaths = getListOfRootFilePaths(path, g_verbose);
    g_outputPath = "/data/local1/raedler/J-PET/Gate_Multi-detector_Post-processing/cmake-build-default/Output/";

    g_mapCenter = {0., 0., 0.};  // mm
    g_mapHalfSize = {0.5, 0.5, 1270.};  // mm
    g_nVoxels = {1, 1, 2540};  // 1×1×1 mm spacing

    runSequentially(g_fullPaths.size(), processSingleFile);
    //runInSeparateProcesses(g_fullPaths.size(), processSingleFile, 128);


    // make the sensitivity map
    // put in the draft
    // add the other event selections
    // make the castor translation too

    return 0;
}