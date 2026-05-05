#include <iostream>
#include <TSystem.h>
#include <fstream>
#include <set>
#include "include/utils.h"
#include "include/eventSelection.h"
#include "include/CASToRTools.h"


// Globals
bool g_verbose;
TString g_treeName;
std::vector<TString> g_fullPaths;
TString g_gantry;
TString g_selection;
TString g_lutName;
std::vector<LutEntry> g_lut;
TString g_outputPath;
TString g_outputFileName;
uint32_t g_timeInitial_ms, g_timeFinal_ms;
// TString g_dzMax_mm;
double g_dzMax_mm;
float g_min_dt_ps, g_max_dt_ps;
TString g_tof_fwhm_ps;


void processSingleFile(const size_t idx) {
    std::cout << "Processing: " << g_fullPaths[idx] << std::endl;

    TFile* file = getTFile(g_fullPaths[idx], "UPDATE", g_verbose);
    TTree* tree = getTTree(file, g_treeName, g_verbose);

    TString fileName = gSystem->BaseName(g_fullPaths[idx]);
    fileName.Remove(fileName.Last('.'));

    if (!checkIfLeafExists(tree, "castorID")) {
        std::cerr << "Error: castorID leaf(s) missing." << std::endl;
        std::exit(1);
    }

    if (!checkIfLeafExists(tree, "trueness")) {
        std::cerr << "Error: trueness leaf missing." << std::endl;
        std::exit(1);
    }

    if (!checkIfLeafExists(tree, "group")) {
        std::cerr << "Error: coincidence grouping leafs missing." << std::endl;
        std::exit(1);
    }

    // Add boolean branch for the preselection
    Bool_t selection;
    TBranch* b = tree->Branch("selection", &selection, "selection/O");
    if (g_selection == "energy") {selectBasedOnEnergy(tree, selection, b, g_verbose);}
    if (g_selection == "time") {selectBasedOnTime(tree, selection, b, g_verbose);}

    // // Add boolean branch for the preselection
    // Bool_t withinMaxAxialDifference;
    // TBranch* b2 = tree->Branch("withinMaxAxialDifference", &withinMaxAxialDifference, "withinMaxAxialDifference/O");
    // checkMaxAxialDifference(tree, withinMaxAxialDifference, b2, g_lut, g_dzMax_mm, g_verbose);

    TLeaf* time1 = tree->GetLeaf("time1");
    TLeaf* time2 = tree->GetLeaf("time2");

    TLeaf* energy1 = tree->GetLeaf("energy1");
    TLeaf* energy2 = tree->GetLeaf("energy2");

    TLeaf* castorID1 = tree->GetLeaf("castorID1");
    TLeaf* castorID2 = tree->GetLeaf("castorID2");

    TLeaf* gantryID1 = tree->GetLeaf("gantryID1");
    TLeaf* gantryID2 = tree->GetLeaf("gantryID2");

    TLeaf* trueness = tree->GetLeaf("trueness");

    std::vector<CdfEntry> buffer;

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t ii = 0; ii < nEntries; ++ii) {
        tree->GetEntry(ii);

        // Both need to be above 200 keV
        if ((energy1->GetValue() < .2) || (energy2->GetValue() < .2)) {continue;}

        if (g_selection == "true") {
            if (!trueness->GetValue()) {continue;}
        } else if (g_selection == "energy") {
            if (!selection) {continue;}
        } else if (g_selection == "time") {
            if (!selection) {continue;}
        }

        // if (!withinMaxAxialDifference) {continue;}
        //if (withinMaxAxialDifference) {continue;}

        TString currentEntryGantryName =  assignGantryName(gantryID1->GetValue(), gantryID2->GetValue());
        if (currentEntryGantryName != g_gantry && g_gantry != "ALL") {continue;}

        double t1 = time1->GetValue();
        double t2 = time2->GetValue();

        uint32_t time_ms = static_cast<uint32_t>(std::round(1000.0 * (t1 + t2) / 2.0));
        float delta_time_ps = static_cast<float>(std::round((t1 - t2) / 1e-12));
        uint32_t c1 = castorID1->GetValue();
        uint32_t c2 = castorID2->GetValue();

        buffer.push_back({time_ms, delta_time_ps, c1, c2});
    }

    // std::cout << buffer.size() << std::endl;

    if (gSystem->AccessPathName(g_outputPath)) {gSystem->mkdir(g_outputPath, true);}
    std::ofstream out(g_outputPath + fileName + ".cdf" , std::ios::binary);
    out.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(CdfEntry));
    out.close();
}



void getGlobalMinMaxTime(std::vector<CdfEntry> castorData) {
    // Smallest value
    auto min_it = std::min_element(castorData.begin(), castorData.end(),
        [](const CdfEntry &a, const CdfEntry &b) { return a.time_ms < b.time_ms; });

    // Largest value
    auto max_it = std::max_element(castorData.begin(), castorData.end(),
        [](const CdfEntry &a, const CdfEntry &b) { return a.time_ms < b.time_ms; });

    if (min_it != castorData.end()) {
        if (min_it->time_ms < g_timeInitial_ms) {
            g_timeInitial_ms = min_it->time_ms;
        }
    }
    if (max_it != castorData.end()) {
        if (max_it->time_ms > g_timeFinal_ms) {
            g_timeFinal_ms = max_it->time_ms;
        }
    }
}



void getGlobalMinMaxDT(std::vector<CdfEntry> castorData) {
    // Smallest value
    auto min_it = std::min_element(castorData.begin(), castorData.end(),
        [](const CdfEntry &a, const CdfEntry &b) { return a.delta_time_ps < b.delta_time_ps; });

    // Largest value
    auto max_it = std::max_element(castorData.begin(), castorData.end(),
        [](const CdfEntry &a, const CdfEntry &b) { return a.delta_time_ps < b.delta_time_ps; });

    if (min_it != castorData.end()) {
        if (min_it->delta_time_ps < g_min_dt_ps) {
            g_min_dt_ps = min_it->delta_time_ps;
        }
    }
    if (max_it != castorData.end()) {
        if (max_it->delta_time_ps > g_max_dt_ps) {
            g_max_dt_ps = max_it->delta_time_ps;
        }
    }
}



size_t mergeCastorDataFiles() {
    size_t nEntriesTotal = 0;

    std::ofstream out(g_outputPath + g_outputFileName + ".cdf" , std::ios::binary);
    for (int ii = 0; ii < g_fullPaths.size(); ++ii) {
        TString fileName = gSystem->BaseName(g_fullPaths[ii]);
        fileName.Remove(fileName.Last('.'));
        auto buffer = readCdfFile(g_outputPath + fileName + ".cdf");

        getGlobalMinMaxTime(buffer);
        getGlobalMinMaxDT(buffer);

        out.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(CdfEntry));
        nEntriesTotal += buffer.size();
        gSystem->Unlink(g_outputPath +  fileName + ".cdf");
    }
    out.close();

    return nEntriesTotal;
}



void makeCastorDataHeaderFile(size_t nEntriesTotal) {
    std::ofstream file(g_outputPath + g_outputFileName + ".cdh" , std::ios::binary);
    if (!file) {throw std::runtime_error("Cannot open header file.");}

    file << "Scanner name: " << g_lutName << "\n";
    file << "Data filename: " << g_outputFileName << ".cdf\n";
    file << "Data type: PET\n";
    file << "Data mode: list-mode\n";  // histogram, normalization
    file << "Isotope: unknown\n";
    file << "Number of events: " << nEntriesTotal << "\n";
    file << "Start time (s): " << g_timeInitial_ms / 1000 << "\n";
    file << "Duration (s): " << (g_timeFinal_ms - g_timeInitial_ms) / 1000 << "\n";
    file << "Maximum number of lines per event: 1\n";
    file << "Calibration factor: 1\n";
    file << "Attenuation correction flag: 0\n";
    file << "Normalization correction flag: 0\n";
    file << "Scatter correction flag: 0\n";
    file << "Random correction flag: 0\n";
    file << "TOF information flag: 1\n";
    file << "TOF resolution (ps): " << g_tof_fwhm_ps << "\n";
    file << "Per event TOF resolution flag: 0\n";
    file << "List TOF measurement range (ps): " << g_max_dt_ps - g_min_dt_ps << "\n";
    // file << "Maximum axial difference mm: " << "-1, -1, -1, -1" << "\n";  // [mm] only for sensitivity image computation for list-mode data and when no normalization file is provided
    // file << "Histo TOF number of bins: " << 10 << "\n";  // required if TOF information flag is on for histogram data
    // file << "Histo TOF bin size (ps): " << 10 << "\n";  // required if TOF information flag is on for histogram data
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

    g_gantry = "ALL";
    g_selection = "true";
    g_lutName =  "TB_J-PET_7th_gen_brain_insert_WHR_4_18_1_mm";
    TString dzMax_mm = "1000";
    g_tof_fwhm_ps = "400";

    std::map<std::string, ArgumentOptions> argOpts = {
        {"gantry", {{"ALL", "TB-TB", "TB-BI", "BI-BI"}, g_gantry}},
        {"selection", {{"true", "energy", "time"}, g_selection}},
        {"lut", {{"TB_J-PET_7th_gen_brain_insert_WHR_4_18_1_mm", "TB_J-PET_7th_gen_brain_insert_WHR_6_30_1_mm"}, g_lutName}},
        {"dzMax", {{"0", "2000"}, dzMax_mm}},
        {"tof", {{"0", "1000"}, g_tof_fwhm_ps}}
    };

    parseArguments(argc, argv, argOpts);
    g_dzMax_mm = std::stod(dzMax_mm.Data());
    checkArguments(argOpts);

    g_outputPath = "/data/local1/raedler/J-PET/Gate_Multi-detector_Post-processing/cmake-build-default/Output/";  // needs to have trailing "/"
    g_outputFileName = g_gantry + "_" + g_selection;
    // g_outputFileName = g_gantry + "_" + g_selection + "_not_" + dzMax_mm;
    g_timeInitial_ms = std::numeric_limits<uint32_t>::max();
    g_timeFinal_ms = std::numeric_limits<uint32_t>::min();
    g_min_dt_ps = std::numeric_limits<float>::max();
    g_max_dt_ps = -std::numeric_limits<float>::max();
    g_lut = readLutBinary("/data/local1/raedler/J-PET/CASToR/castor/config/scanner/" + g_lutName + ".lut");

    // Adapt the same folder structure (last two folders) from the input
    std::vector<TString> pathSplit = splitPath(path);
    g_outputPath += pathSplit[pathSplit.size() - 2] + "/" + pathSplit[pathSplit.size() - 1] + "/";

    // runSequentially(g_fullPaths.size(), processSingleFile);
    runInSeparateProcesses(g_fullPaths.size(), processSingleFile, 128);

    size_t nEntriesTotal = mergeCastorDataFiles();
    // size_t nEntriesTotal = 1;

    makeCastorDataHeaderFile(nEntriesTotal);

    // make the sensitivity map
    // put in the draft

    return 0;
}