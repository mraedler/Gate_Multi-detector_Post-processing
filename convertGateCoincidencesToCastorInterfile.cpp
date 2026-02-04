#include <iostream>
#include <TSystem.h>
#include <fstream>
#include "include/utils.h"
#include "include/CASToRTools.h"


// Globals
bool g_verbose;
TString g_treeName;
std::vector<TString> g_fullPaths;
TString g_lut_name;
TString g_outputPath;
TString g_outputFileName;
uint32_t g_timeInitial_ms, g_timeFinal_ms;
float g_min_dt_ps, g_max_dt_ps;
TString g_tof_fwhm_ps;


void processSingleFile(const size_t idx) {
    std::cout << "Processing: " << g_fullPaths[idx] << std::endl;

    TFile* file = getTFile(g_fullPaths[idx], "READ", g_verbose);
    TTree* tree = getTTree(file, g_treeName, g_verbose);

    TString fileName = gSystem->BaseName(g_fullPaths[idx]);
    fileName.Remove(fileName.Last('.'));

    if (!checkIfLeafExists(tree, "castorID1") || !checkIfLeafExists(tree, "castorID2")) {
        std::cerr << "Error: castorID leaf(s) missing." << std::endl;
        std::exit(1);
    }

    if (!checkIfLeafExists(tree, "trueness")) {
        std::cerr << "Error: trueness leaf missing." << std::endl;
        std::exit(1);
    }

    TLeaf* time1 = tree->GetLeaf("time1");
    TLeaf* time2 = tree->GetLeaf("time2");

    TLeaf* castorID1 = tree->GetLeaf("castorID1");
    TLeaf* castorID2 = tree->GetLeaf("castorID2");

    TLeaf* trueness = tree->GetLeaf("trueness");

    std::vector<CdfEntry> buffer;

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t ii = 0; ii < nEntries; ++ii) {
        tree->GetEntry(ii);

        // True events only
        if (!trueness->GetValue()) {continue;}

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

    std::exit(1);
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

    file << "Scanner name: " << g_lut_name << "\n";
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
    g_lut_name = "TB_J-PET_7th_gen_brain_insert_dz_1_mm";
    g_outputPath = "/data/local1/raedler/J-PET/Gate_Multi-detector_Post-processing/cmake-build-default/Output/";  // needs to have trailing "/"
    g_outputFileName = "true";
    g_timeInitial_ms = std::numeric_limits<uint32_t>::max();
    g_timeFinal_ms = std::numeric_limits<uint32_t>::min();
    g_min_dt_ps = std::numeric_limits<float>::max();
    g_max_dt_ps = -std::numeric_limits<float>::max();
    g_tof_fwhm_ps = "500";

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
    // add the other event selections

    return 0;
}