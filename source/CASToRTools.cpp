#include "../include/CASToRTools.h"
#include "../include/utils.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <TF1.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TGraph.h>
#include <TApplication.h>
#include <TCanvas.h>

double g_FWHM_sigma_conversion = 2.0 * std::sqrt(2.0 * std::log(2.0));

std::vector<LutEntry> readLutBinary(const char* filePath)
{
    std::ifstream file(filePath, std::ios::binary | std::ios::ate);
    if (!file) {
        throw std::runtime_error("Cannot open LUT file!");
    }

    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    if (size % sizeof(LutEntry) != 0) {
        throw std::runtime_error("LUT file size mismatch!");
    }

    std::vector<LutEntry> lut(size / sizeof(LutEntry));
    file.read(reinterpret_cast<char*>(lut.data()), size);

    return lut;
}



ScannerParams totalBodyJPETWithBrainInsert_4_18()
{
	return {
		        {
			        {2, 2, 24, 16, 330},
					{2, 3, 24, 16, 600},
					{2, 1, 12, 16, 330}
		        },
				{
		            {0.0, 6.0 / g_FWHM_sigma_conversion},
					{0.0, 6.0 / g_FWHM_sigma_conversion},
					{0.0, 4.0 / g_FWHM_sigma_conversion}
				},
				1.0
			};
}



ScannerParams totalBodyJPETWithBrainInsert_6_30()
{
	return {
		        {
			        {2, 2, 24, 16, 330},
					{2, 3, 24, 16, 600},
					{2, 1, 12, 11, 330}
		        },
				{
		            {0.0, 6.0 / g_FWHM_sigma_conversion},
					{0.0, 6.0 / g_FWHM_sigma_conversion},
					{0.0, 6.0 / g_FWHM_sigma_conversion}
				},
				1.0
			};
}



ScannerParams GEDiscoveryMI()
{
	return {
			        {
				        {34, 4, 4, 9, 4}
			        },
					{
			            {0.0, 0.0 / g_FWHM_sigma_conversion}
					},
					0.0
				};
}



void setCastorID(TTree* tree, const TLeaf* gantryID, const TLeaf* rsectorID, const TLeaf* crystalID, const TLeaf* layerID, Int_t& castorID, TBranch* b, TString& lutName) {
	// Get the scanner parameters
	ScannerParams sp;
	if (lutName == "TB_J-PET_7th_gen_brain_insert_WHR_6_30_1_mm") {
		sp = totalBodyJPETWithBrainInsert_6_30();
	} else if (lutName == "TB_J-PET_7th_gen_brain_insert_WHR_4_18_1_mm") {
		sp = totalBodyJPETWithBrainInsert_4_18();
	} else {
		std::cerr << "Unknown LUT name." << std::endl;
		std::exit(1);
	}

	// Set the distributions and the random number generator for the blurring along the axial direction
	std::vector<std::normal_distribution<>> dists;
	dists.reserve(sp.gaussParams.size());
	for (const auto& v : sp.gaussParams) {
		dists.emplace_back(v[0], v[1]);
	}
	std::mt19937 gen(std::random_device{}());

	// Gantry dependent shift
	std::vector<std::uint32_t> shift = {0};
	std::uint32_t cumulative = 0;
	for (const auto& row : sp.gantryShape) {
		std::uint32_t prod = 1;
		for (auto v : row) {
			prod *= v;
		}
		cumulative += prod;
		shift.push_back(cumulative);
	}

	// Should print: [      0  506880 1889280 2016000] for the totalBodyJPETWithBrainInsert_4_18 parameters
	//printVectorEntries(shift);

	std::vector<std::uint32_t> layerEntries;
	for (const auto& row : sp.gantryShape) {
		layerEntries.push_back(row[3] * row[4]);
	}

	// Should print: [5280 9600 5280] for the totalBodyJPETWithBrainInsert_4_18 parameters
	//printVectorEntries(layerEntries);

	//
	Long64_t nEntries = tree->GetEntries();
	// Long64_t nEntries = tree->GetEntries() / 100;

	// Allocate the output vector
	// std::vector<int> castorIDs(nEntries);

	for (Long64_t ii = 0; ii < nEntries; ++ii) {
		tree->GetEntry(ii);
		int gID = gantryID->GetValue();

		// Extract the layer from the layerID
		int layerNumber = layerID->GetValue() / layerEntries[gID]; // Integer division, i.e. including floor
		int layerIdx = layerID->GetValue() - layerNumber * layerEntries[gID];

		// Disentangle the moduleID and the longitudinalID from the layerIdx
		// Unravel multi-index in column-major (Fortran) order
		int moduleID = layerIdx % sp.gantryShape[gID][3];
		int longitudinalID = layerIdx / sp.gantryShape[gID][3];

		// Add gaussian noise along the longitudinal dimension
		int counter = 0;
		int randomShift;
		bool exceedsLimits;
		do {
			counter++;
			randomShift = static_cast<int>(std::round(dists[gID](gen) / sp.longitudinalSpacing));
			exceedsLimits = longitudinalID + randomShift < 0 || longitudinalID + randomShift >= sp.gantryShape[gID][4];
			//if (exceedsLimits) std::cout << "Trial " << counter << " rejected\n";
		} while (exceedsLimits);
		longitudinalID += randomShift;

		// Fortran-style index linearlization (undoing the above separation)
		layerIdx = longitudinalID * sp.gantryShape[gID][3] + moduleID;

		// Ravel multi-index in row-major C-style in four dimensions
		// with indices (i0, i1, i2, i3) = (layerNumber, crystalID, rsectorID, layerIdx)
		// and shape (d0, d1, d2, d3) depending on the gantry
		// flat_idx = ((i0 ∗ d1 + i1) ∗ d2 + i2) ∗ d3 + i3
		castorID = ((layerNumber * sp.gantryShape[gID][1] + crystalID->GetValue()) * sp.gantryShape[gID][2] + rsectorID->GetValue()) * layerEntries[gID] + layerIdx + shift[gID];
		b->Fill();
	}
}



void setCastorID(TTree* tree, const TLeaf* rsectorID, const TLeaf* moduleID, const TLeaf* submoduleID, const TLeaf* crystalID, const TLeaf* globalPosZ, Int_t& castorID, TBranch* b, TString& lutName) {
	// Get the scanner parameters
	ScannerParams sp;
	if (lutName == "GE_Discovery_MI") {
		sp = GEDiscoveryMI();
	} else {
		std::cerr << "Unknown LUT name." << std::endl;
		std::exit(1);
	}
	//
	Long64_t nEntries = tree->GetEntries();

	TGraph *gr = new TGraph();
	std::mt19937 gen(std::random_device{}());
	std::normal_distribution<float> normal(0., 0.1);

	for (Long64_t ii = 0; ii < nEntries; ++ii) {
		tree->GetEntry(ii);
		// {34, 4, 4, 9, 4}
		int rsID = rsectorID->GetValue();
		int mID = moduleID->GetValue();
		int smID = submoduleID->GetValue();
		int cID = crystalID->GetValue();

		// Un-ravel the crystal ID (C-style)
		int cID_zed = cID / sp.gantryShape[0][4];
		int cID_lat = cID % sp.gantryShape[0][4];

		// Ravel for consistency check
		// int cID2 = cID_zed * sp.gantryShape[0][4] + cID_lat;
		// std::cout << cID - cID2  << std::endl;

		if (mID == 0) {
			gr->SetPoint(gr->GetN(), cID_zed + normal(gen), globalPosZ->GetValue());
		}
		// C-style raveling
		castorID = ((rsID * sp.gantryShape[0][1] + mID)
		                  * sp.gantryShape[0][2] + smID)
		                  * sp.gantryShape[0][3] * sp.gantryShape[0][4] + cID;
		b->Fill();
	}
	delete gr;
	// TApplication app("app", 0, nullptr);
	// TCanvas canvas("c", "c", 1200, 600);
	// gr->SetMarkerStyle(24); // filled circle
	// gr->Draw("AP");
	// app.Run();
}



void checkCastorID(TTree* tree, const TLeaf* gantryID, const TLeaf* globalPosX, const TLeaf* globalPosY, const TLeaf* globalPosZ, Int_t& castorID, const std::vector<LutEntry>& lut) {
	// const ScannerParams sp = totalBodyJPETWithBrainInsert_4_18();
	const ScannerParams sp = GEDiscoveryMI();
	int nGantries = sp.gantryShape.size();

	//
	std::vector<std::vector<float>> depthDeviation(nGantries);
	std::vector<std::vector<float>> lateralDeviation(nGantries);
	std::vector<std::vector<float>> longitudinalDeviation(nGantries);

	Long64_t nEntries = tree->GetEntries();
	auto nEntriesFloat = static_cast<float>(nEntries);
	for (Long64_t ii = 0; ii < nEntries; ++ii) {
		tree->GetEntry(ii);

		int gID = gantryID->GetValue();

		// Rotate
		float v_x = lut[castorID].OrVx;
		float v_y = lut[castorID].OrVy;
		float lut_depth   = lut[castorID].Posx * v_x + lut[castorID].Posy * v_y;
		float lut_lateral = lut[castorID].Posx * v_y - lut[castorID].Posy * v_x;
		float depth   = globalPosX->GetValue() * v_x + globalPosY->GetValue() * v_y;
		float lateral = globalPosX->GetValue() * v_y - globalPosY->GetValue() * v_x;

		depthDeviation[gID].push_back(depth - lut_depth);
		lateralDeviation[gID].push_back(lateral - lut_lateral);
		longitudinalDeviation[gID].push_back(globalPosZ->GetValue() - lut[castorID].Posz);
	}

	// Print the min/max deviation
	std::cout << std::endl;
	std::cout << std::fixed << std::setprecision(1) << "Min, max deviation [mm]\n=======================" << std::endl;
	for (int kk = 0; kk < nGantries; ++kk) {
		std::cout << "Gantry " << kk << "\n--------" << std::endl;
		auto [minDepthDeviation, maxDepthDeviation] = minmax(depthDeviation[kk]);
		std::cout << "Depth:        " << minDepthDeviation << ", " << maxDepthDeviation << std::endl;

		auto [minLateralDeviation, maxLateralDeviation] = minmax(lateralDeviation[kk]);
		std::cout << "Lateral:      " << minLateralDeviation << ", " << maxLateralDeviation << std::endl;

		auto [minLongitudinalDeviation, maxLongitudinalDeviation] = minmax(longitudinalDeviation[kk]);
		std::cout << "Longitudinal: " << minLongitudinalDeviation << ", " << maxLongitudinalDeviation << std::endl;

		std::cout << std::endl;
	}

	// Visualize and check consistency
	TApplication app("app", 0, nullptr);
	TCanvas canvas("c", "c", 1200, 600);
	canvas.Divide(3, nGantries);

	std::cout << std::fixed << std::setprecision(1) << "Counts within histogram\n=======================" << std::endl;

	for (int jj = 0; jj < nGantries; jj++) {

		canvas.cd(1 + jj * nGantries);
		// TH1D* histDepthDeviation = getHistogram(depthDeviation[jj], -16., 1., 32, std::to_string(1 + jj * nGantries).c_str(), "; Depth deviation; Count");
		TH1D* histDepthDeviation = getHistogram(depthDeviation[jj], -16., 1., 32, std::to_string(1 + jj * nGantries).c_str(), "; Depth deviation; Count");
		float itg = histDepthDeviation->Integral();
		std::cout << "Depth:        " << static_cast<int>(itg) << " / " << nEntries << " (" << itg / nEntriesFloat * 100 << " %)" << std::endl;
		histDepthDeviation->Draw();
		histDepthDeviation->SetStats(0);

		// TF1* f = new TF1("f", "expo", -15, 15);
		// f->SetParameters(11., -1.5e-02);
		// f->SetLineColor(kRed);
		// f->SetLineWidth(2);
		// f->Draw("SAME");

		// histDepthDeviation->Fit(f);

		canvas.cd(2 + jj * nGantries);
		// TH1D* histLateralDeviation = getHistogram(lateralDeviation[jj], -4., .5, 16, std::to_string(2 + jj * nGantries).c_str(), "; Lateral deviation; Count");
		TH1D* histLateralDeviation = getHistogram(lateralDeviation[jj], -8., 1., 16, std::to_string(2 + jj * nGantries).c_str(), "; Lateral deviation; Count");
		itg = histLateralDeviation->Integral();
		std::cout << "Lateral:      " << static_cast<int>(itg) << " / " << nEntries << " (" << itg / nEntriesFloat * 100 << " %)" << std::endl;
		histLateralDeviation->Draw();
		histLateralDeviation->SetStats(0);

		canvas.cd(3 + jj * nGantries);
		TH1D* histLongitudinalDeviation = getHistogram(longitudinalDeviation[jj], -8., 1., 16, std::to_string(3 + jj * nGantries).c_str(), "; Longitudinal deviation; Count");
		itg = histLongitudinalDeviation->Integral();
		std::cout << "Longitudinal: " << static_cast<int>(itg) << " / " << nEntries << " (" << itg / nEntriesFloat * 100 << " %)" << std::endl;
		histLongitudinalDeviation->Draw();
		histLongitudinalDeviation->SetStats(0);
	}
	app.Run();

}



std::vector<CdfEntry> readCdfFile(const TString& filename) {
	std::ifstream in(filename, std::ios::binary);
	if (!in) {throw std::runtime_error("Cannot open file: " + filename);}

	// Determine file size
	in.seekg(0, std::ios::end);
	std::streamsize fileSize = in.tellg();
	in.seekg(0, std::ios::beg);

	if (fileSize % sizeof(CdfEntry) != 0) {throw std::runtime_error("File size is not a multiple of CdfEntry size");}

	size_t nEntries = fileSize / sizeof(CdfEntry);

	std::vector<CdfEntry> data(nEntries);

	if (!in.read(reinterpret_cast<char*>(data.data()), fileSize)) {throw std::runtime_error("Error reading file: " + filename);}

	return data;
}


std::vector<CdfEntryWithoutTOF> readCdfWithoutTOFFile(const TString& filename) {
	std::ifstream in(filename, std::ios::binary);
	if (!in) {throw std::runtime_error("Cannot open file: " + filename);}

	// Determine file size
	in.seekg(0, std::ios::end);
	std::streamsize fileSize = in.tellg();
	in.seekg(0, std::ios::beg);

	if (fileSize % sizeof(CdfEntryWithoutTOF) != 0) {throw std::runtime_error("File size is not a multiple of CdfEntry size");}

	size_t nEntries = fileSize / sizeof(CdfEntryWithoutTOF);

	std::vector<CdfEntryWithoutTOF> data(nEntries);

	if (!in.read(reinterpret_cast<char*>(data.data()), fileSize)) {throw std::runtime_error("Error reading file: " + filename);}

	return data;
}