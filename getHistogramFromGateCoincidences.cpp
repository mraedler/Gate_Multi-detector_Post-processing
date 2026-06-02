#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TH3D.h>
#include <THn.h>
#include <TSystem.h>
#include <TRandom3.h>
#include "TROOT.h"
#include "include/vec3.h"
#include "include/utils.h"
#include "include/eventSelection.h"

// Globals
bool g_verbose;
TString g_treeName;
std::vector<TString> g_fullPaths;
TString g_outputPath;

int g_nDims = 0;
std::vector<int> g_nVoxels;
std::vector<double> g_mapCenter;
std::vector<double> g_mapHalfSize;
std::vector<double> g_xMin;
std::vector<double> g_xMax;



double getPositronRangeSample(TRandom3& rng) {
    constexpr double c = 0.516;
    constexpr double k_1 = 37.9;  // [1/mm]
    constexpr double k_2 = 3.10;  // [1/mm]

    constexpr double tau_1 = 1. / k_1;
    constexpr double tau_2 = 1. / k_2;

    constexpr double w_1 = c / k_1 / (c / k_1 + (1 - c) / k_2);
    // constexpr double w_2 = (1 - c) / k_2 / (c / k_1 + (1 - c) / k_2);  // w_2 = 1 - w_1

    double sign = (rng.Uniform() < 0.5) ? -1. : 1.;

    if (rng.Uniform() < w_1) {
        return rng.Exp(tau_1) * sign;
    } else {
        return rng.Exp(tau_2) * sign;
    }
}


void processSingleFile(const size_t idx) {
    std::cout << "Processing: " << g_fullPaths[idx] << std::endl;

    TFile* file = getTFile(g_fullPaths[idx], "UPDATE", g_verbose);
    TTree* tree = getTTree(file, g_treeName, g_verbose);
    // listAvailableBranches(tree);

    TString fileName = gSystem->BaseName(g_fullPaths[idx]);
    fileName.Remove(fileName.Last('.'));

    if (gSystem->AccessPathName(g_outputPath)) {gSystem->mkdir(g_outputPath, true);}

    // Allocate the output histograms

    // TH3D* h_TBTB = new TH3D(fileName + "_TBTB", "Counts",
    //                    g_nVoxels[0], g_mapCenter[0] - g_mapHalfSize[0], g_mapCenter[0] + g_mapHalfSize[0],
    //                    g_nVoxels[1], g_mapCenter[1] - g_mapHalfSize[1], g_mapCenter[1] + g_mapHalfSize[1],
    //                    g_nVoxels[2], g_mapCenter[2] - g_mapHalfSize[2], g_mapCenter[2] + g_mapHalfSize[2]);

    // TH3D* h_TBBI = new TH3D(fileName + "_TBBI", "Counts",
    //                    g_nVoxels[0], g_mapCenter[0] - g_mapHalfSize[0], g_mapCenter[0] + g_mapHalfSize[0],
    //                    g_nVoxels[1], g_mapCenter[1] - g_mapHalfSize[1], g_mapCenter[1] + g_mapHalfSize[1],
    //                    g_nVoxels[2], g_mapCenter[2] - g_mapHalfSize[2], g_mapCenter[2] + g_mapHalfSize[2]);

    // TH3D* h_BIBI = new TH3D(fileName + "_BIBI", "Counts",
    //                    g_nVoxels[0], g_mapCenter[0] - g_mapHalfSize[0], g_mapCenter[0] + g_mapHalfSize[0],
    //                    g_nVoxels[1], g_mapCenter[1] - g_mapHalfSize[1], g_mapCenter[1] + g_mapHalfSize[1],
    //                    g_nVoxels[2], g_mapCenter[2] - g_mapHalfSize[2], g_mapCenter[2] + g_mapHalfSize[2]);

    THnD* h_TBTB = new THnD(fileName + "_TBTB", "Counts", g_nDims, g_nVoxels.data(), g_xMin.data(), g_xMax.data());
    THnD* h_TBBI = new THnD(fileName + "_TBBI", "Counts", g_nDims, g_nVoxels.data(), g_xMin.data(), g_xMax.data());
    THnD* h_BIBI = new THnD(fileName + "_BIBI", "Counts", g_nDims, g_nVoxels.data(), g_xMin.data(), g_xMax.data());

    // std::exit(1);

    // Add boolean branch for the preselection
    Bool_t selection;
    TBranch* b = tree->Branch("selection", &selection, "selection/O");
    // selectBasedOnTime(tree, selection, b, g_verbose);
    selectBasedOnEnergy(tree, selection, b, g_verbose);

    TLeaf* trueness = tree->GetLeaf("trueness");

    TLeaf* energy1 = tree->GetLeaf("energy1");
    TLeaf* energy2 = tree->GetLeaf("energy2");

    TLeaf* gantryID1 = tree->GetLeaf("gantryID1");
    TLeaf* gantryID2 = tree->GetLeaf("gantryID2");

    TLeaf* crystalID1 = tree->GetLeaf("crystalID1");
    TLeaf* crystalID2 = tree->GetLeaf("crystalID2");

    TLeaf* sourcePosX1 = tree->GetLeaf("sourcePosX1");
    TLeaf* sourcePosY1 = tree->GetLeaf("sourcePosY1");
    TLeaf* sourcePosZ1 = tree->GetLeaf("sourcePosZ1");
    // TLeaf* sourcePosX2 = tree->GetLeaf("sourcePosX2");
    // TLeaf* sourcePosY2 = tree->GetLeaf("sourcePosY2");
    // TLeaf* sourcePosZ2 = tree->GetLeaf("sourcePosZ2");

    TLeaf* globalPosX1 = tree->GetLeaf("globalPosX1");
    TLeaf* globalPosY1 = tree->GetLeaf("globalPosY1");
    TLeaf* globalPosZ1 = tree->GetLeaf("globalPosZ1");

    TLeaf* globalPosX2 = tree->GetLeaf("globalPosX2");
    TLeaf* globalPosY2 = tree->GetLeaf("globalPosY2");
    TLeaf* globalPosZ2 = tree->GetLeaf("globalPosZ2");

    int nDims = g_nVoxels.size();

    TRandom3 rng(12345);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t ii = 0; ii < nEntries; ++ii) {
        tree->GetEntry(ii);
        bool aboveEnergyThreshold = (energy1->GetValue() > .2) && (energy2->GetValue() > .2);  // both above 200 keV

        // if (trueness->GetValue() && aboveEnergyThreshold) {
        // if (selection && aboveEnergyThreshold) {
        if (selection && trueness->GetValue() && aboveEnergyThreshold) {

            // if (!((gantryID1->GetValue() == 0) & (gantryID2->GetValue() == 0) & ((crystalID1->GetValue() == 0) | (crystalID1->GetValue() == 1)) & (crystalID2->GetValue() == 0))) continue;

            //
            Vec3 A{globalPosX1->GetValue(), globalPosY1->GetValue(), globalPosZ1->GetValue()};
            Vec3 B{globalPosX2->GetValue(), globalPosY2->GetValue(), globalPosZ2->GetValue()};
            // Vec3 S{sourcePosX1->GetValue(), sourcePosY1->GetValue(), sourcePosZ1->GetValue()};
            Vec3 S{sourcePosX1->GetValue() + getPositronRangeSample(rng),
                sourcePosY1->GetValue() + getPositronRangeSample(rng),
                sourcePosZ1->GetValue() + getPositronRangeSample(rng)};

            double lor_length = norm(B - A);
            // std::cout << lor_length << std::endl;
            double source_lor_distance = norm(cross(S - A, B - A)) / lor_length;
            // if (source_lor_distance > 1e-4) continue;

            // double t = dot(S - A, B - A) / norm2(B - A);
            double t = dot(S - A, B - A) / (lor_length * lor_length);
            double p_tilde = 2. * (t - 0.5);
            double alpha = 0.5 / 360. * 2. * 3.14159265358979323846;
            double fwhm = lor_length / 2 * alpha / 2 * (1 - p_tilde * p_tilde);

            // std::cout << t << std::endl;

            // todo
            double values[nDims] = {sourcePosX1->GetValue(), sourcePosY1->GetValue(), sourcePosZ1->GetValue()};
            // double values[nDims] = {source_lor_distance};
            // double values[nDims] = {getPositronRangeSample(rng)};
            // double values[nDims] = {t};
            // double values[nDims] = {lor_length};
            // double values[nDims] = {fwhm};
            // double values[nDims] = {fwhm * fwhm};
            // double values[nDims] = {1 / fwhm};
            // double values[nDims] = {lor_length, p_tilde};

            TString currentEntryGantryName = assignGantryName(gantryID1->GetValue(), gantryID2->GetValue());

            if (currentEntryGantryName == "TB-TB") {
                // h_TBTB->Fill(sourcePosX1->GetValue(), sourcePosY1->GetValue(), sourcePosZ1->GetValue());
                h_TBTB->Fill(values);
            } else if (currentEntryGantryName == "TB-BI") {
                // h_TBBI->Fill(sourcePosX1->GetValue(), sourcePosY1->GetValue(), sourcePosZ1->GetValue());
                h_TBBI->Fill(values);
            } else if (currentEntryGantryName == "BI-BI") {
                // h_BIBI->Fill(sourcePosX1->GetValue(), sourcePosY1->GetValue(), sourcePosZ1->GetValue());
                h_BIBI->Fill(values);
            } else {
                std::cerr << "Warning: unknown gantry.\n";
            }
        }
    }

    h_TBTB->SaveAs(g_outputPath + fileName + "_TBTB.root");
    h_TBBI->SaveAs(g_outputPath + fileName + "_TBBI.root");
    h_BIBI->SaveAs(g_outputPath + fileName + "_BIBI.root");

    delete h_TBTB;
    delete h_TBBI;
    delete h_BIBI;

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



THnD* getTHnD(const TString& fileName) {
    TFile* file = TFile::Open(g_outputPath +  fileName + ".root", "READ");
    THnD* hist = (THnD*)file->Get(fileName);
    THnD* histCopy = (THnD*)hist->Clone();
    file->Close();
    delete file;
    gSystem->Unlink(g_outputPath +  fileName + ".root");
    return histCopy;
}



void mergeHists(const TString& type) {
    TString firstFileName = gSystem->BaseName(g_fullPaths[0]);
    firstFileName.Remove(firstFileName.Last('.'));

    // TH3D* firstHist = getTH3D(firstFileName + "_" + type);
    THnD* firstHist = getTHnD(firstFileName + "_" + type);

    for (unsigned int ii = 1; ii < g_fullPaths.size(); ii++) {
        TString fileName = gSystem->BaseName(g_fullPaths[ii]);
        fileName.Remove(fileName.Last('.'));

        // TH3D* hist = getTH3D(fileName + "_" + type);
        THnD* hist = getTHnD(fileName + "_" + type);
        firstHist->Add(hist);
        delete hist;
    }
    // firstHist->SetName("TH3");
    firstHist->SetName("THnD");
    // firstHist->SaveAs(g_outputPath + type + "_true.root");
    // firstHist->SaveAs(g_outputPath + type + "_p_distribution.root");
    // firstHist->SaveAs(g_outputPath + type + "_2L_distribution.root");
    // firstHist->SaveAs(g_outputPath + type + "_2L_p_distribution.root");
    // firstHist->SaveAs(g_outputPath + type + "_FWHM_distribution.root");
    // firstHist->SaveAs(g_outputPath + type + "_FWHM_squared_distribution.root");
    // firstHist->SaveAs(g_outputPath + type + "_FWHM_inverse_distribution.root");
    // firstHist->SaveAs(g_outputPath + type + "_time.root");
    // firstHist->SaveAs(g_outputPath + type + "_energy.root");
    // firstHist->SaveAs(g_outputPath + type + "_time_true.root");
    firstHist->SaveAs(g_outputPath + type + "_energy_true.root");
    // firstHist->SaveAs(g_outputPath + type + "_d_min.root");
    // firstHist->SaveAs(g_outputPath + type + "_d_min_plus_range.root");
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
    g_outputPath = "/data/local1/raedler/J-PET/Gate_Multi-detector_Post-processing/cmake-build-default/Output/";  // needs to have trailing "/"
    std::vector<TString> pathSplit = splitPath(path);
    g_outputPath += pathSplit[pathSplit.size() - 2] + "/" + pathSplit[pathSplit.size() - 1] + "/";

    // 1×1×1 mm spacing
    g_nVoxels = {1, 1, 2540};
    g_mapCenter = {0., 0., 0.};  // mm
    g_mapHalfSize = {0.5, 0.5, 1270.};  // mm

    // g_nVoxels = {100};
    // g_mapCenter = {0.};
    // g_mapHalfSize = {1e-4};

    // // p-distribution
    // g_nVoxels = {101};
    // g_mapCenter = {0.5};
    // g_mapHalfSize = {0.5};

    // // 2L-distribution
    // g_nVoxels = {200};
    // g_mapCenter = {1000};
    // g_mapHalfSize = {1000};

    // // fwhm-distribution
    // g_nVoxels = {400};
    // g_mapCenter = {2.};
    // g_mapHalfSize = {2.};

    // // fwhm-squared-distribution
    // g_nVoxels = {400};
    // g_mapCenter = {1.};
    // g_mapHalfSize = {1.};

    // // 2L-p-distribution
    // g_nVoxels = {200, 201};
    // g_mapCenter = {1000., 0.};
    // g_mapHalfSize = {1000., 1.};

    // //
    // g_nVoxels = {100};
    // g_mapCenter = {2.5};
    // g_mapHalfSize = {2.5};


    g_nDims = g_nVoxels.size();

    //
    for (int ii = 0; ii < g_nDims; ii++) {
        g_xMin.push_back(g_mapCenter[ii] - g_mapHalfSize[ii]);
        g_xMax.push_back(g_mapCenter[ii] + g_mapHalfSize[ii]);
    }

    // runSequentially(g_fullPaths.size(), processSingleFile);
    runInSeparateProcesses(g_fullPaths.size(), processSingleFile, 128);

    mergeHists("TBTB");
    mergeHists("TBBI");
    mergeHists("BIBI");

    return 0;
}