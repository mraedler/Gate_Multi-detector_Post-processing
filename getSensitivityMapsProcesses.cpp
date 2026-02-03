#include "ROOT/TThreadExecutor.hxx"
#include "ROOT/RDataFrame.hxx"
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TSystem.h>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TH3.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLine.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <random>
#include <typeinfo>

#include <unistd.h>   // fork(), _exit()
#include <sys/wait.h> // wait()
#include "include/utils.h"



void printTH3BinEdges(const TH3& h) {
	int nx = h.GetNbinsX();
	int ny = h.GetNbinsY();
	int nz = h.GetNbinsZ();

	/*
	std::cout << "X bin edges: ";
	for (int ii = 1; ii <= nx+1; ii++)  // ROOT bins are 1-based
		std::cout << h.GetXaxis()->GetBinLowEdge(ii) << " ";
	std::cout << "\n";

	std::cout << "Y bin edges: ";
	for (int i = 1; i <= ny+1; ++i)
		std::cout << h.GetYaxis()->GetBinLowEdge(i) << " ";
	std::cout << "\n";
	*/

	std::cout << "Z bin edges: ";
	for (int i = 1; i <= nz+1; ++i)
		std::cout << h.GetZaxis()->GetBinLowEdge(i) << " ";
	std::cout << "\n";
}


/*
struct LutEntry {
    float Posx;
    float Posy;
    float Posz;
    float OrVx;
    float OrVy;
    float OrVz;
};
*/

/*
std::vector<LutEntry> readLutBinary(const char* filePath)
{
    std::ifstream file(filePath, std::ios::binary | std::ios::ate);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<LutEntry> lut(size / sizeof(LutEntry));

    file.read(reinterpret_cast<char*>(lut.data()), size);

    return lut;
}
*/


template <typename T>
void printVectorEntries(std::vector<T>& arr) {
	for (int ii = 0; ii < arr.size(); ii++) {
		std::cout << arr[ii] << std::endl;
	}
	std::cout << std::endl;
}


/*
template <typename T>
TH1D* getHistogram(std::vector<T>& arr, const float lowerEdge, const float binWidth, const int nBins, const char* name, const char* labeling) {

	TH1D* hist = new TH1D(name, labeling, nBins, lowerEdge, lowerEdge + nBins * binWidth);

    for (const T& a : arr)
        hist->Fill(a);

	return hist;
}
*/

/*
std::vector<int> getCASToRID(TTree* tree, TLeaf* gantryID, TLeaf* rsectorID, TLeaf* crystalID, TLeaf* layerID, TLeaf* globalPosX, TLeaf* globalPosY, TLeaf* globalPosZ, const std::vector<LutEntry>& lut, bool vis) {
	constexpr std::array<std::array<std::uint32_t, 5>, 3> gantryShape{{
		{2, 2, 24, 16, 330},
		{2, 3, 24, 16, 600},
		{2, 1, 12, 16, 330}}};
	
	int nGantries = gantryShape.size();
	
	// Random sampling for the z blurring; the FWHMs are in [mm]
	std::mt19937 gen(std::random_device{}());
	double FWHM_sigma_conversion = 2.0 * std::sqrt(2.0 * std::log(2.0));
	std::vector<std::normal_distribution<double>> dists = {
		std::normal_distribution<double>{0.0, 6.0 / FWHM_sigma_conversion},
		std::normal_distribution<double>{0.0, 6.0 / FWHM_sigma_conversion},
		std::normal_distribution<double>{0.0, 4.0 / FWHM_sigma_conversion}
	};
	double dz = 1.;  // [mm]

	//
	std::vector<std::uint32_t> shift = {0};
	std::uint32_t cumulative = 0;
	for (const auto& row : gantryShape) {
		std::uint32_t prod = 1;
		for (auto v : row) {
			prod *= v;
		}
		cumulative += prod;
		shift.push_back(cumulative);
	}
	
	// Should print: [      0  506880 1889280 2016000]
	//printVectorEntries(shift);
	
	std::vector<std::uint32_t> layerEntries;
	for (const auto& row : gantryShape) {
		layerEntries.push_back(row[3] * row[4]);
	}
	
	// Should print: [5280 9600 5280]
	//printVectorEntries(layerEntries);
	
	// 
	//Long64_t nEntries = tree->GetEntries() / 100;
	Long64_t nEntries = tree->GetEntries();
	int gID, layerNumber, layerIdx, moduleID, longitudinalID, counter, randomShift, CASToRID;
	float v_x, v_y, lut_depth, lut_lateral, depth, lateral;
	bool exceedsLimits;
	
	// For the vis
	std::vector<std::vector<float>> depthDeviation(nGantries);
	std::vector<std::vector<float>> lateralDeviation(nGantries);
	std::vector<std::vector<float>> longitudinalDeviation(nGantries);
	
	// Allocate the output vector
	std::vector<int> CASToRIDs(nEntries);
	
	for (Long64_t ii = 0; ii < nEntries; ++ii) {
		tree->GetEntry(ii);
		gID = gantryID->GetValue();
		
		// Extract the layer from the layerID
		layerNumber = layerID->GetValue() / layerEntries[gID]; // Integer division, i.e. including floor
		layerIdx = layerID->GetValue() - layerNumber * layerEntries[gID];
		
		// Disentangle the moduleID and the longitudinalID from the layerIdx
		// Unravel multi-index in column-major (Fortran) order
		moduleID = layerIdx % gantryShape[gID][3];
		longitudinalID = layerIdx / gantryShape[gID][3];
		
		// Add gaussian noise along the longitidunal dimension
		counter = 0;
		do {
			counter++;
			randomShift = static_cast<int>(std::round(dists[gID](gen) / dz));
			exceedsLimits = longitudinalID + randomShift < 0 || longitudinalID + randomShift >= gantryShape[gID][4];
			//if (exceedsLimits) std::cout << "Trial " << counter << " rejected\n";
		} while (exceedsLimits);
		longitudinalID += randomShift;
		
		// Fortran-style index linearlization (undoing the above separation)
		layerIdx = longitudinalID * gantryShape[gID][3] + moduleID;
		
		// Ravel multi-index in row-major C-style in four dimensions 
		// with indices (i0, i1, i2, i3) = (layerNumber, crystalID, rsectorID, layerIdx)
		// and shape (d0, d1, d2, d3) depending on the gantry
		// flat_idx = ((i0 ∗ d1 + i1) ∗ d2 + i2) ∗ d3 + i3
		CASToRID = ((layerNumber * gantryShape[gID][1] + crystalID->GetValue()) * gantryShape[gID][2] + rsectorID->GetValue()) * layerEntries[gID] + layerIdx + shift[gID];
		CASToRIDs[ii] = CASToRID;
		
		// 
		if (vis) {
			v_x = lut[CASToRID].OrVx;
			v_y = lut[CASToRID].OrVy;
			lut_depth   = lut[CASToRID].Posx * v_x + lut[CASToRID].Posy * v_y;
			lut_lateral = lut[CASToRID].Posx * v_y - lut[CASToRID].Posy * v_x;
			depth   = globalPosX->GetValue() * v_x + globalPosY->GetValue() * v_y;
			lateral = globalPosX->GetValue() * v_y - globalPosY->GetValue() * v_x;
			
			depthDeviation[gID].push_back(depth - lut_depth);
			lateralDeviation[gID].push_back(lateral - lut_lateral);
			longitudinalDeviation[gID].push_back(globalPosZ->GetValue() - lut[CASToRID].Posz);
		}
		
	}
	
	if (vis) {
		// Visualize and check consistency
		TApplication app("app", 0, nullptr);
		TCanvas canvas("c", "c", 1200, 600);
		canvas.Divide(3, nGantries);
		
		for (int jj = 0; jj < nGantries; jj++) {
			
			canvas.cd(1 + jj * nGantries);
			TH1D* histDepthDeviation = getHistogram(depthDeviation[jj], -16., 1., 32, std::to_string(1 + jj * nGantries).c_str(), "; Depth deviation; Count");
			histDepthDeviation->Draw();
			histDepthDeviation->SetStats(0);
			
			canvas.cd(2 + jj * nGantries);
			TH1D* histLateralDeviation = getHistogram(lateralDeviation[jj], -4., .5, 16, std::to_string(2 + jj * nGantries).c_str(), "; Lateral deviation; Count");
			histLateralDeviation->Draw();
			histLateralDeviation->SetStats(0);
			
			canvas.cd(3 + jj * nGantries);
			TH1D* histLongitudinalDeviation = getHistogram(longitudinalDeviation[jj], -8., .1, 160, std::to_string(3 + jj * nGantries).c_str(), "; Longitudinal deviation; Count");
			histLongitudinalDeviation->Draw();
			histLongitudinalDeviation->SetStats(0);
		}
		app.Run();
	}
    
    return CASToRIDs;
}
*/

/*
std::vector<bool> setMinSectorDifference(TTree* tree, const TLeaf* gantryID1, const TLeaf* gantryID2, const TLeaf* rsectorID1, const TLeaf* rsectorID2, int minSectorDifference, bool plotHistogram) {
	
	//Long64_t nEntries = tree->GetEntries() / 100;
	Long64_t nEntries = tree->GetEntries();
	
	std::vector<int> sectorDifference(nEntries, -1);
	//std::vector<bool> aboveMinSectorDifference(nEntries, true);
	std::vector<bool> aboveMinSectorDifference(nEntries, false);
	
	int gID1, gID2, rsID1, rsID2, absolute_sector_difference;
	float passingCounter = 0.;
	
	for (Long64_t ii = 0; ii < nEntries; ii++) {
		tree->GetEntry(ii);
		
		gID1 = gantryID1->GetValue();
		gID2 = gantryID2->GetValue();
		
		rsID1 = rsectorID1->GetValue();
		rsID2 = rsectorID2->GetValue();
		
		absolute_sector_difference = std::abs(rsID2 - rsID1);

		if (gID1 == 0 and gID2 == 0) {
			sectorDifference[ii] = std::min(absolute_sector_difference, 24 - absolute_sector_difference);
		} else if (gID1 == 1 and gID2 == 1) {
			sectorDifference[ii] = std::min(absolute_sector_difference, 24 - absolute_sector_difference);
		} else if (gID1 == 2 and gID2 == 2) {
			sectorDifference[ii] = std::min(absolute_sector_difference, 12 - absolute_sector_difference);
		}
		
		//if (sectorDifference[ii] >= 0 and sectorDifference[ii] < minSectorDifference) {
		//	aboveMinSectorDifference[ii] = false;
		//}
		if (sectorDifference[ii] == -1 or sectorDifference[ii] >= minSectorDifference) {
			aboveMinSectorDifference[ii] = true;
			passingCounter++;
		}
	}
	
	if (plotHistogram) {
		std::cout << "Passing percentage minimum sector difference: " << passingCounter / nEntries * 100 << " %." << std::endl;
	
		TApplication app("app", 0, nullptr);
		TH1D* histSectorDifference = getHistogram(sectorDifference, -2.5, 1, 15, "h", "; Sector difference; Count");
		TCanvas canvas("c", "c", 800, 600);
		histSectorDifference->Draw();
		histSectorDifference->SetStats(0);
		app.Run();
		
	}
	
	return aboveMinSectorDifference;
}
*/

/*
double findFirstMinimumAfterZero(std::vector<double> scatterTestForHistogram, const bool vis) {
	// 
	float lowerEdgeHist = -100.5;  // [cm]
	float binWidthHist = 1.;  // [cm]
	int nBins = 400 + 1;
	
	int idx_0 = static_cast<int>(std::round(-(lowerEdgeHist + binWidthHist / 2) / binWidthHist + 1));
	
	TH1D* histScatterTest = getHistogram(scatterTestForHistogram, lowerEdgeHist, binWidthHist, 401, "h", "; Scatter test [cm]; Count");
	
	int thresholdIndex;
	
	for (int jj = idx_0; jj <= histScatterTest->GetNbinsX() - 1; ++jj) {
		//double binsCenter = (histScatterTest->GetBinCenter(jj + 1) + histScatterTest->GetBinCenter(jj)) / 2;
		double countDifference = histScatterTest->GetBinContent(jj + 1) - histScatterTest->GetBinContent(jj);
		//std::cout << jj << "  binsCenter = " << binsCenter << "  countDifference = " << countDifference << std::endl;
		
		// Search for the point where the derivative changes
		if (countDifference > 0.) {
			thresholdIndex = jj;
			break;
		}
		
	}
	
	double binCenterMinimum = histScatterTest->GetBinCenter(thresholdIndex);
	//std::cout << binCenterMinimum << std::endl;
	
	if (vis) {
		TApplication app("app", 0, nullptr);
						
		TCanvas canvas("c", "c", 800, 600);
		histScatterTest->Draw();
		histScatterTest->SetStats(0);
		
		canvas.Update();
		TLine* line = new TLine(binCenterMinimum, gPad->GetUymin(), binCenterMinimum, gPad->GetUymax());
		//line->SetLineColor(kRed);
		//line->SetLineWidth(2);
		line->Draw("same");
		
		app.Run();
	}
	
	return binCenterMinimum;
	
}
*/

/*
std::vector<bool> runScatterTest(TTree* tree, TLeaf* time1, TLeaf* time2, TLeaf* gantryID1, TLeaf* gantryID2, std::vector<LutEntry> lut, std::vector<int> CASToRID1, std::vector<int> CASToRID2, std::vector<bool> aboveMinSectorDifference, bool verbose) {
	
	//Long64_t nEntries = tree->GetEntries() / 100;
	Long64_t nEntries = tree->GetEntries();
	
	float x1, y1, z1, x2, y2, z2;
	double t1, t2;
	
	double distance;
	double speedOfLight = 2.99792458e11;  // [mm / s]
	
	std::vector<double> scatterTest(nEntries);
	std::vector<double> scatterTestForHistogram;
	
	for (Long64_t ii = 0; ii < nEntries; ++ii) {
		
		tree->GetEntry(ii);
		
		x1 = lut[CASToRID1[ii]].Posx;
		y1 = lut[CASToRID1[ii]].Posy;
		z1 = lut[CASToRID1[ii]].Posz;
		
		x2 = lut[CASToRID2[ii]].Posx;
		y2 = lut[CASToRID2[ii]].Posy;
		z2 = lut[CASToRID2[ii]].Posz;
		
		t1 = time1->GetValue();
		t2 = time2->GetValue();
		
		distance = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));  // [mm]
		scatterTest[ii] = (distance - (t2 - t1) * speedOfLight) / 10.;  // [cm]
		
		if (aboveMinSectorDifference[ii]) {
			scatterTestForHistogram.push_back(scatterTest[ii]);
		}
	}
	
	double threshold = findFirstMinimumAfterZero(scatterTestForHistogram, verbose);
	//std::cout << threshold << std::endl;
	
	// different thresholds for different gantries
	double thresholdTBTB = 25.;  // [cm]
	double thresholdOther = 10.;  // [cm]
	
	std::vector<bool> passingScatterTest(nEntries, false);
	float passingCounter = 0.;
	
	for (Long64_t ii = 0; ii < nEntries; ++ii) {
		tree->GetEntry(ii);
		
		if ((gantryID1->GetValue() < 2) && (gantryID2->GetValue() < 2)) {
			if (scatterTest[ii] >= thresholdTBTB) {
				passingScatterTest[ii] = true;
				passingCounter++;
			}
		} else {
			if (scatterTest[ii] >= thresholdOther) {
				passingScatterTest[ii] = true;
				passingCounter++;
			}
		}
	}
	
	if (verbose) {std::cout << "Passing percentage scatter test: " << passingCounter / nEntries * 100 << " %." << std::endl;}
	
	return passingScatterTest;
}
*/


std::vector<bool> andVector(std::vector<bool>& a, std::vector<bool>& b, std::string abc, bool verbose) {
	
	if (a.size() !=  b.size()) {
		throw std::runtime_error("Vectors a and b must have the same size");
	}
	
	Long64_t n = a.size();
	
	std::vector<bool> c(n);
	
	float trueCounter = 0;
	
	for (Long64_t ii = 0; ii < n; ii++) {
		c[ii] = a[ii] && b[ii];
		trueCounter += c[ii];
	}
	
	if (verbose) {std::cout << "Passing percentage " << abc << ": "  << trueCounter / n * 100 << " %." << std::endl;}
	//if (verbose) {std::cout << "Passing percentage " << abc << ": "  << trueCounter << " %." << std::endl;}
	
	return c;
};



template <typename T>
std::vector<T> getSelectionVector(TTree* tree, TLeaf* leaf, const std::vector<bool>& selection) {
	std::vector<T> result;
    //result.reserve(selection.size());

    for (Long64_t ii = 0; ii < selection.size(); ++ii) {
        if (!selection[ii]) continue;

        tree->GetEntry(ii);
        result.push_back(static_cast<T>(leaf->GetValue()));
    }

    return result;
}


/*
std::tuple<std::vector<Long64_t>, std::vector<Long64_t>, std::vector<Long64_t>, std::vector<Long64_t>> groupCoincidences(TTree* tree, TLeaf* time1, TLeaf* time2, TLeaf* energy1, TLeaf* energy2, std::vector<bool> preSelection)
{
	//Long64_t nEntries = tree->GetEntries() / 100;
	Long64_t nEntries = tree->GetEntries();

	Long64_t ii = 0;

	std::vector<Long64_t> groupEdges = {0};
	std::vector<Long64_t> groupMultiplicities;
	std::vector<Long64_t> idx1, idx2;

	while (ii < nEntries) {

		//std::set<std::tuple<double, double>> seen;
		std::map<std::tuple<double, double>, int> seen;
		int nextId = 0;

		while (ii < nEntries) {
			tree->GetEntry(ii);

			auto k1 = std::make_tuple(time1->GetValue(), energy1->GetValue());
			auto k2 = std::make_tuple(time2->GetValue(), energy2->GetValue());

			// If seen is not empty and both times have not been seen
			if (!seen.empty() && !seen.count(k1) && !seen.count(k2)) {
				break;
			}

			// The std::set only inserts values (or tuples here) that it has not seen before; so no additional check necessary here
			//seen.insert(k1);
			//seen.insert(k2);

			// Using the std::map instead
			if (!seen.count(k1)) {
				seen[k1] = nextId++;
			}
			if (!seen.count(k2)) {
				seen[k2] = nextId++;
			}

			idx1.push_back(seen[k1]);
			idx2.push_back(seen[k2]);

			ii++;
		}

		groupEdges.push_back(ii);
		groupMultiplicities.push_back(seen.size());
	}

	if (groupEdges.back() != nEntries) {
		std::cerr << "Warning: not all elements grouped!\n";
	}


	// Print all time entries
	for (Long64_t jj = 0; jj < nEntries; jj++) {
		tree->GetEntry(jj);
		std::cout << std::setprecision(16) << jj << " " << time1->GetValue() << " " << time2->GetValue() << std::endl;
	}
	std::cout << std::endl;



	// Print grouped time entries
	for (Long64_t jj = groupEdges.size() - 100; jj < groupEdges.size() - 1; jj++){
		std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;
		for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; kk++) {
			tree->GetEntry(kk);
			std::cout << std::setprecision(16) << time1->GetValue() << " " << time2->GetValue() << std::endl;
		}
		std::cout << std::endl;
	}



	std::cout << groupEdges.size() << std::endl;
	std::cout << groupMultiplicities.size() << std::endl;
	std::cout << idx1.size() << std::endl;
	std::cout << idx2.size() << std::endl;
	std::cout << nEntries << std::endl;



	// Check for multiplicity 1, indicating whether the criterion to identfy identical events is insufficient
	for (Long64_t jj = 0; jj < groupEdges.size() - 1; jj++){
		if (groupMultiplicities[jj] == 1) {
			std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;
			for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; kk++) {
				tree->GetEntry(kk);
				//std::cout << std::setprecision(16) << time1->GetValue() << " " << time2->GetValue() << std::endl;
				std::cout << std::setprecision(16) << time1->GetValue() - time2->GetValue() << std::endl;
				std::cout << std::setprecision(16) << energy1->GetValue() - energy2->GetValue() << std::endl;
			}
			std::cout << std::endl;
		}
	}



	// Check the indexing
	for (Long64_t jj = groupEdges.size() - 100; jj < groupEdges.size() - 1; jj++){
		std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;
		for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; kk++) {
			tree->GetEntry(kk);
			std::cout << std::setprecision(16) << time1->GetValue() << " " << idx1[kk] << " " << time2->GetValue() << " " << idx2[kk] << std::endl;
		}
		std::cout << std::endl;
	}


	return {groupEdges, groupMultiplicities, idx1, idx2};
}
*/

/*
std::tuple<std::vector<Long64_t>, std::vector<Long64_t>, std::vector<Long64_t>, std::vector<Long64_t>> groupCoincidences2(std::vector<double>& time1, std::vector<double>& time2, std::vector<float>& energy1, std::vector<float>& energy2)
{	
	Long64_t nEntries = time1.size();
	Long64_t ii = 0;
	
	std::vector<Long64_t> groupEdges = {0};
	std::vector<Long64_t> groupMultiplicities;
	std::vector<Long64_t> idx1, idx2;

	while (ii < nEntries) {

		//std::set<std::tuple<double, double>> seen;
		std::map<std::tuple<double, float>, int> seen;
		int nextId = 0;

		while (ii < nEntries) {
			
			auto k1 = std::make_tuple(time1[ii], energy1[ii]);
			auto k2 = std::make_tuple(time2[ii], energy2[ii]);
			
			// If seen is not empty and both times have not been seen
			if (!seen.empty() && !seen.count(k1) && !seen.count(k2)) {
				break;
			}
			
			// The std::set only inserts values (or tuples here) that it has not seen before; so no additional check necessary here
			//seen.insert(k1);
			//seen.insert(k2);
			
			// Using the std::map instead
			if (!seen.count(k1)) {
				seen[k1] = nextId++;
			}
			if (!seen.count(k2)) {
				seen[k2] = nextId++;
			}

			idx1.push_back(seen[k1]);
			idx2.push_back(seen[k2]);
			
			ii++;
		}
		
		groupEdges.push_back(ii);
		groupMultiplicities.push_back(seen.size());
	}
	
	if (groupEdges.back() != nEntries) {
		std::cerr << "Warning: not all elements grouped!\n";
	}
	

	// Print all time entries
	for (Long64_t jj = 0; jj < nEntries; jj++) {
		std::cout << std::setprecision(16) << jj << " " << time1[jj] << " " << time2[jj] << std::endl;
	}
	std::cout << std::endl;

	

	// Print grouped time entries
	for (Long64_t jj = groupEdges.size() - 100; jj < groupEdges.size() - 1; jj++){
		std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;
		for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; kk++) {
			std::cout << std::setprecision(16) << time1[kk] << " " << time2[kk] << std::endl;
		}
		std::cout << std::endl;
	}

	

	std::cout << groupEdges.size() << std::endl;
	std::cout << groupMultiplicities.size() << std::endl;
	std::cout << idx1.size() << std::endl;
	std::cout << idx2.size() << std::endl;
	std::cout << nEntries << std::endl;

	

	// Check for multiplicity 1, indicating whether the criterion to identfy identical events is insufficient
	for (Long64_t jj = 0; jj < groupEdges.size() - 1; jj++){
		if (groupMultiplicities[jj] == 1) {
			std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;
			for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; kk++) {
				//std::cout << std::setprecision(16) << time1[kk] << " " << time2[kk] << std::endl;
				std::cout << std::setprecision(16) << time1[kk] - time2[kk] << std::endl;
				std::cout << std::setprecision(16) << energy1[kk] - energy2[kk] << std::endl;
			}
			std::cout << std::endl;
		}
	}

	

	// Check the indexing
	for (Long64_t jj = groupEdges.size() - 100; jj < groupEdges.size() - 1; jj++){
		std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;
		for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; kk++) {
			std::cout << std::setprecision(16) << time1[kk] << " " << idx1[kk] << " " << time2[kk] << " " << idx2[kk] << std::endl;
		}
		std::cout << std::endl;
	}

	
	return {groupEdges, groupMultiplicities, idx1, idx2};
}
*/

/*
void plotGroupStatistics(const std::vector<Long64_t>& groupEdges, std::vector<Long64_t>& groupMultiplicities) {
	// Get the group sizes
	std::vector<Long64_t> groupSizes;
    for (size_t ii = 0; ii < groupEdges.size() - 1; ++ii)
        groupSizes.push_back(groupEdges[ii + 1] - groupEdges[ii]);
	
	// Expected to peak at the triangle numbers: 1, 3, 6, 10, 15, 21, ...
	TApplication app("app", 0, nullptr);
	
	TH1D* histGroupSizes = getHistogram(groupSizes, -0.5, 1, 29, "h1", ";Group size; Count");
	TH1D* histGroupMultiplicities = getHistogram(groupMultiplicities, 0.5, 1, 12, "h2", "; Multiplicity; Count");

    TCanvas canvas("c", "c", 1200, 600);
	canvas.Divide(2,1);
	canvas.cd(1);
    gPad->SetLogy();
    histGroupSizes->Draw();
    histGroupSizes->SetStats(0);
    
	canvas.cd(2);
    gPad->SetLogy();
    histGroupMultiplicities->Draw();
    histGroupMultiplicities->SetStats(0);
	
    app.Run();
}
*/

/*
std::vector<bool> selectBasedOnTime(const std::vector<Long64_t>& groupEdges, const std::vector<Long64_t>& idx1, const std::vector<Long64_t>& idx2, const bool verbose) {
	std::vector<bool> selection(idx1.size(), false);
	
	for (Long64_t jj = 0; jj < groupEdges.size() - 1; jj++){
		if (verbose) {std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;}
		
		std::set<Long64_t> usedEvents;
		for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; kk++) {
			if (verbose) {std::cout << idx1[kk] << " " << idx2[kk] << " ";}
			
			//if (usedEvents.count(idx1[kk]) || usedEvents.count(idx2[kk])) {
			if (usedEvents.count(idx1[kk]) || usedEvents.count(idx2[kk]) || usedEvents.size() > 0) {  // only choose one coincidence
				if (verbose) {std::cout << std::endl;}
				continue;
			}
			
			if (verbose) {std::cout << "*" << std::endl;}
			
			selection[kk] = true;
			usedEvents.insert(idx1[kk]);
			usedEvents.insert(idx2[kk]);
			
		}
		if (verbose) {std::cout << std::endl;}
	}
	
	return selection;
}
*/

/*
std::vector<bool> selectBasedOnEnergy(const std::vector<Long64_t>& groupEdges, const std::vector<Long64_t>& idx1, const std::vector<Long64_t>& idx2, const std::vector<float>& energy1, const std::vector<float>& energy2, const bool verbose) {
	std::vector<bool> selection(idx1.size(), false);
	
	for (Long64_t jj = 0; jj < groupEdges.size() - 1; jj++){
		if (verbose) {std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;}
		
		// Sort the group indices in decending order with respect to the total energy
		std::vector<Long64_t> sortedGroupIndices(groupEdges[jj + 1] - groupEdges[jj]);
		std::iota(sortedGroupIndices.begin(), sortedGroupIndices.end(), groupEdges[jj]);
		std::sort(sortedGroupIndices.begin(), sortedGroupIndices.end(), [&](Long64_t a, Long64_t b) {return (energy1[a] + energy2[a]) > (energy1[b] + energy2[b]);});
		
		// Selection
		std::set<Long64_t> usedEvents;
		for (Long64_t ll = 0; ll < sortedGroupIndices.size(); ll++) {
			int sortedGroupIndex = sortedGroupIndices[ll];
			
			if (verbose) {std::cout << idx1[sortedGroupIndex] << " " << idx2[sortedGroupIndex] << " " << energy1[sortedGroupIndex] + energy2[sortedGroupIndex] << " ";}
			
			//if (usedEvents.count(idx1[sortedGroupIndex]) || usedEvents.count(idx2[sortedGroupIndex])) {
			if (usedEvents.count(idx1[sortedGroupIndex]) || usedEvents.count(idx2[sortedGroupIndex]) || usedEvents.size() > 0) {  // only choose one coincidence
				if (verbose) {std::cout << std::endl;}
				continue;
			}
			
			if (verbose) {std::cout << "*" << std::endl;}
			
			selection[sortedGroupIndex] = true;
			usedEvents.insert(idx1[sortedGroupIndex]);
			usedEvents.insert(idx2[sortedGroupIndex]);
			
		}
		if (verbose) {std::cout << std::endl;}
	
	}
	
	return selection;
}
*/


std::vector<bool> resetSelectionVector(const std::vector<bool>& selectionDense, const std::vector<bool>& selectionSparse) {
	std::vector<bool> selectionDenseSparse(selectionSparse.size(), false);
	
	Long64_t jj = 0;
	
	for (Long64_t ii = 0; ii < selectionSparse.size(); ii++) {
		if (!selectionSparse[ii]) {continue;}
		selectionDenseSparse[ii] = selectionDense[jj];
		jj++;
	}
	
	if (jj != selectionDense.size()) {
		throw std::runtime_error("Number of true elements in the sparse vector does not match the number of elements in the dense vector.");
	}
	
	return selectionDenseSparse;
}



TH3D* getSensitivityMap(const TString& fullPath, const TString& treeName, const TString& gantryName, const std::vector<double>& mapCenter, const std::vector<double>& mapHalfSize, const std::vector<int>& nVoxels)
{
	// Allocate the output array
	TH3D* h = new TH3D(gSystem->BaseName(fullPath), "Counts", 
					   nVoxels[0], mapCenter[0] - mapHalfSize[0], mapCenter[0] + mapHalfSize[0], 
					   nVoxels[1], mapCenter[1] - mapHalfSize[1], mapCenter[1] + mapHalfSize[1],
					   nVoxels[2], mapCenter[2] - mapHalfSize[2], mapCenter[2] + mapHalfSize[2]);
	//return h;

	//printTH3BinEdges(*h);

	TFile f(fullPath, "READ");
	if (f.IsZombie()) {
		std::cerr << "Error opening file " << fullPath << std::endl;
		return h;
	}
	listAvailableTrees(&f);

	TTree* tree = (TTree*)f.Get(treeName);
	if (!tree) {
		std::cerr << "Tree " << treeName << " not found in " << fullPath << std::endl;
		return h;
    }

	listAvailableLeaves(tree);

	// For the loop over all entries
	//Long64_t nEntries = tree->GetEntries() / 100;
	Long64_t nEntries = tree->GetEntries();
	std::cout << gSystem->BaseName(fullPath) << ": " << nEntries << " entries" << std::endl;

	//std::cout << nVoxels[0] << std::endl;

	// Get relevant leafs
	TLeaf* comptonCrystal1 = tree->GetLeaf("comptonCrystal1");
	TLeaf* comptonCrystal2 = tree->GetLeaf("comptonCrystal2");
	TLeaf* rayleighCrystal1 = tree->GetLeaf("RayleighCrystal1");
	TLeaf* rayleighCrystal2 = tree->GetLeaf("RayleighCrystal2");
	TLeaf* eventID1 = tree->GetLeaf("eventID1");
	TLeaf* eventID2 = tree->GetLeaf("eventID2");
	TLeaf* time1 = tree->GetLeaf("time1");
	TLeaf* time2 = tree->GetLeaf("time2");
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
	
	TLeaf* rsectorID1 = tree->GetLeaf("rsectorID1");
	TLeaf* rsectorID2 = tree->GetLeaf("rsectorID2");
	
	TLeaf* crystalID1 = tree->GetLeaf("crystalID1");
	TLeaf* crystalID2 = tree->GetLeaf("crystalID2");
	
	TLeaf* layerID1 = tree->GetLeaf("layerID1");
	TLeaf* layerID2 = tree->GetLeaf("layerID2");
	
	TLeaf* globalPosX1 = tree->GetLeaf("globalPosX1");
	TLeaf* globalPosX2 = tree->GetLeaf("globalPosX2");
	
	TLeaf* globalPosY1 = tree->GetLeaf("globalPosY1");
	TLeaf* globalPosY2 = tree->GetLeaf("globalPosY2");
	
	TLeaf* globalPosZ1 = tree->GetLeaf("globalPosZ1");
	TLeaf* globalPosZ2 = tree->GetLeaf("globalPosZ2");
	
	auto lut = readLutBinary("/data/local1/raedler/J-PET/CASToR/castor/config/scanner/TB_J-PET_7th_gen_brain_insert_dz_1_mm.lut");
	std::vector<int> CASToRID1 = getCASToRID(tree, gantryID1, rsectorID1, crystalID1, layerID1, globalPosX1, globalPosY1, globalPosZ1, lut, false);
	std::vector<int> CASToRID2 = getCASToRID(tree, gantryID2, rsectorID2, crystalID2, layerID2, globalPosX2, globalPosY2, globalPosZ2, lut, false);
	
	std::vector<bool> aboveMinSectorDifference = setMinSectorDifference(tree, gantryID1, gantryID2, rsectorID1, rsectorID2, 2, false);
	std::vector<bool> passingScatterTest = runScatterTest(tree, time1, time2, gantryID1, gantryID2, lut, CASToRID1, CASToRID2, aboveMinSectorDifference, false);
	std::vector<bool> preSelection = andVector(aboveMinSectorDifference, passingScatterTest, "pre-selection", false);
	
	// Only pass the pre-selection to the grouping stage
	auto time1S = getSelectionVector<double>(tree, time1, preSelection);
	auto time2S = getSelectionVector<double>(tree, time2, preSelection);
	auto energy1S = getSelectionVector<float>(tree, energy1, preSelection);
	auto energy2S = getSelectionVector<float>(tree, energy2, preSelection);
	
	// For the event selection: group coincidences that use the same single
	//auto [groupEdges, groupMultiplicities, idx1, idx2] = groupCoincidences(tree, time1, time2, energy1, energy2, preSelection);
	auto [groupEdges, groupMultiplicities, idx1, idx2] = groupCoincidences2(time1S, time2S, energy1S, energy2S);
	//plotGroupStatistics(groupEdges, groupMultiplicities);
	
	// Select with removing duplicate entries
	// Based on time (first) or based on energy or based on truth
	
	//auto selectionTime = selectBasedOnTime(groupEdges, idx1, idx2, false);
	auto selectionEnergy = selectBasedOnEnergy(groupEdges, idx1, idx2, energy1S, energy2S, false);
	
	//std::vector<bool> selectionTimeReset = resetSelectionVector(selectionTime, preSelection);
	std::vector<bool> selectionEnergyReset = resetSelectionVector(selectionEnergy, preSelection);
	
	//std::vector<bool> completeSelection = andVector(preSelection, selectionTimeReset, "event selection time", false);
	std::vector<bool> completeSelection = andVector(preSelection, selectionEnergyReset, "event selection energy", false);
	
	//float counttt = 0;
	//float cc2 = 0;
	//float gg = 0;
	//float abc = 0;

	bool sameEvent, firstCompton, zeroRayleigh, trueCoincidence, aboveThreshold, inGantry;

	//for (Long64_t ii = 0; ii < nEntries / 100; ii++) {
	for (Long64_t ii = 0; ii < nEntries; ii++) {
		tree->GetEntry(ii);
		sameEvent = eventID1->GetValue() == eventID2->GetValue();
		firstCompton = (comptonCrystal1->GetValue() == 1) && (comptonCrystal2->GetValue() == 1);
		zeroRayleigh = (rayleighCrystal1->GetValue() == 0) && (rayleighCrystal2->GetValue() == 0);
		aboveThreshold = (energy1->GetValue() > .2) && (energy2->GetValue() > .2);  // both above 200 keV
		trueCoincidence = sameEvent && firstCompton && zeroRayleigh;
		
		//if (lol[ii] & trueCoincidence) {abc++;}
		//if (ff[ii] & trueCoincidence) {counttt++;}
		//if (trueCoincidence) {cc2++;}
		//if (!preSelection[ii] && trueCoincidence) {gg++;}
		
		
		// With the 2 and 3 ring TB-J-PET geometry, so that there are 3 gantries with the brain insert included
		//std::cout << gantryID1->GetValue() << " " << gantryID2->GetValue() << std::endl;
		if (gantryName == "Comb.") {
			inGantry = true;
		} else if (gantryName == "TB-TB") {
			inGantry = (gantryID1->GetValue() < 2) && (gantryID2->GetValue() < 2);  
		} else if (gantryName == "TB-BI") {
			inGantry = ((gantryID1->GetValue() < 2) && (gantryID2->GetValue() == 2)) || ((gantryID1->GetValue() == 2) && (gantryID2->GetValue() < 2));
		} else if (gantryName == "BI-BI") {
			inGantry = (gantryID1->GetValue() == 2) && (gantryID2->GetValue() == 2);
		} else {
			std::cerr << "Error: unknown gantry.\n";
			std::exit(EXIT_FAILURE);
		}

		//if (trueCoincidence && aboveThreshold && inGantry) {
		//if (completeSelection[ii] && aboveThreshold && inGantry) {
		if (completeSelection[ii] && trueCoincidence && aboveThreshold && inGantry) {
			//std::cout << gSystem->BaseName(fullPath) << " true found" << std::endl;
			h->Fill(sourcePosX1->GetValue(), sourcePosY1->GetValue(), sourcePosZ1->GetValue()); 
		}
	}
	
	//std::cout << counttt << std::endl;
	//std::cout << abc << std::endl;
	//std::cout << cc2 << std::endl;
	//std::cout << gg << std::endl;
	//std::cout << groupEdges.size() - 1 << std::endl;
	
	return h;
}

/*
void getSensitivityMapAndSave(const TString fullPath, const TString& treeName, const std::vector<double>& mapCenter, const std::vector<double>& mapHalfSize, const std::vector<int>& nVoxels)
{
    TH3D* h = getSensitivityMap(fullPath, treeName, mapCenter, mapHalfSize, nVoxels);

    TFile f(outFile, "RECREATE");
    if (h) h->Write();
    f.Close();

    delete h;
}
*/


TH3D* mergeSensitivityMaps(TH3D* a, TH3D* b) {
	a->Add(b);
	delete b;
	return a;
}



auto getSensitivityMapsThreads(const std::vector<TString>& fullPaths, const std::string& treeName, const std::string& gantryName, const std::vector<double>& mapCenter, const std::vector<double>& mapHalfSize, const std::vector<int>& nVoxels)
{
	// Run in parallel
	int nThreads = 16;  // (todo: remove hard coding)
	ROOT::EnableImplicitMT(nThreads);
	ROOT::TThreadExecutor pool;

	// Keep individual
	auto h_maps = pool.Map([&](const TString fullPath) -> TH3D* {
		return getSensitivityMap(fullPath, treeName, gantryName, mapCenter, mapHalfSize, nVoxels);
		}, fullPaths);
	
	/*
	// Merge on the fly
	TH3D* h_map = pool.MapReduce([&](const TString fullPath) -> TH3D* {
		return getSensitivityMap(fullPath, treeName, gantryName, mapCenter, mapHalfSize, nVoxels);
		}, fullPaths, mergeSensitivityMaps);
	*/
	
	return h_maps;
}


/*
void getSensitivityMapsProcesses(const std::vector<TString>& fullPaths, const std::string& treeName, const std::string& gantryName, const std::vector<double>& mapCenter, const std::vector<double>& mapHalfSize, const std::vector<int>& nVoxels, const std::string& outputDir)
{
    int maxProcesses = 128;
    std::vector<pid_t> children;
  

    for (size_t ii = 0; ii < fullPaths.size(); ++ii) {
        pid_t pid = fork();
        //std::cout << pid << std::endl; 
        
        if (pid == 0) {
            // CHILD process
            //TString outFile = Form("%s_output.root", fullPaths[i].Data());
            TH3D* h_map = getSensitivityMap(fullPaths[ii], treeName, gantryName, mapCenter, mapHalfSize, nVoxels);
            h_map->SaveAs((outputDir + "/" + gSystem->BaseName(fullPaths[ii])).c_str());
            _exit(0);
        } else if (pid > 0) {
            // PARENT process
            children.push_back(pid);

            // limit concurrent workers
            if ((int)children.size() >= maxProcesses) {
                int status;
                wait(&status);
                children.pop_back();
            }
        } else {
            std::cerr << "fork() failed!" << std::endl;
            exit(1);
        }
        
    }

    // wait for remaining children
    int status;
    while (wait(&status) > 0);
}
*/


void printTH3D(TH3D* h) {
	int nx = h->GetNbinsX();
	int ny = h->GetNbinsY();
	int nz = h->GetNbinsZ();
    
	//std::cout << nx << " " << ny  << " " << nz << std::endl;

	for (int ix = 1; ix <= nx; ++ix) {
		for (int iy = 1; iy <= ny; ++iy) {
			for (int iz = 1; iz <= nz; ++iz) {
				double content = h->GetBinContent(ix, iy, iz);
				double x = h->GetXaxis()->GetBinCenter(ix);
				double y = h->GetYaxis()->GetBinCenter(iy);
				double z = h->GetZaxis()->GetBinCenter(iz);
                
				//std::cout << "Bin (" << ix << "," << iy << "," << iz << ") "
				//		  << "Center (" << x << "," << y << "," << z << ") "
				//		  << "Content: " << content << std::endl;
            }
        }
    }
}


/*
TH3D* getHist(const TString& fullPath) {
	TFile f(fullPath, "READ");
	std::filesystem::path fullPathFS(fullPath.Data());
	//std::cout << fullPathFS.filename() << std::endl;
	TH3D* h = (TH3D*)f.Get(fullPathFS.filename().c_str());
	h->SetDirectory(0);  // detach from f so that it can be closed and deleted
	f.Close();
	std::filesystem::remove(fullPathFS);
	return h;
}
*/

/*
void mergeAllSensitivityMaps(const std::string& outputDir, const std::string& fileName) {
	// Get list of .root files
	TSystemDirectory dir("dir", outputDir.c_str());
	TList* files = dir.GetListOfFiles();
	if (!files) {
		std::cerr << "Error: cannot open directory " << outputDir << std::endl;
		// return 1;
	}

	// files->Sort();

	// Allocate the array of .root files fullpaths
	std::vector<TString> fullPaths;

	TIter next(files);
	TSystemFile* file;
	while ((file = (TSystemFile*) next())) {
		TString fileName = file->GetName();
			if (!file->IsDirectory() && fileName.BeginsWith("results") && fileName.EndsWith(".root")) {
			fullPaths.push_back(outputDir + '/' + fileName);
		}
	}
	//std::cout << fullPaths.size() << std::endl;
	
	// Get the first before and merge the rest on top.
	TH3D* h_merged = getHist(fullPaths[0]);
	
	for (unsigned int ii = 1; ii < fullPaths.size(); ii++) {
		std::cout << fullPaths[ii] << std::endl;
		
		TH3D* h = getHist(fullPaths[ii]);
		//printTH3D(h);
		
		h_merged->Add(h);
		
	}
	h_merged->SetName("TH3");
	h_merged->SaveAs((outputDir + "/" + fileName + "_true_energy.root").c_str());
}
*/


void processOneOutputDirectory(const std::filesystem::path& outputBaseDir, const std::string& dirPath, const std::string& gantryName, const std::string& treeName, const std::vector<double>& mapCenter, const std::vector<double>& mapHalfSize, const std::vector<int>& nVoxels) {
	// Get list of .root files
	TSystemDirectory dir("dir", dirPath.c_str());
	TList* files = dir.GetListOfFiles();
	
	if (!files) {
		std::cerr << "Error: cannot open directory " << dirPath << std::endl;
	}

	// files->Sort();

	// Allocate the array of .root files fullpaths
	std::vector<TString> fullPaths;

	TIter next(files);
	TSystemFile* file;
	while ((file = (TSystemFile*) next())) {
		TString fileName = file->GetName();
			if (!file->IsDirectory() && fileName.EndsWith(".root")) {
			fullPaths.push_back(dirPath + fileName);
		}
	}
	
	//for (unsigned int ii = 0; ii < fullPaths.size(); ii++) {
	//	std::cout << fullPaths[ii] << std::endl;
	//}
	
	// For the output: get the name of the directory containing the root files
	//std::string basename = getBasename(dirPath);
	std::filesystem::path dirPathFS(dirPath);
	std::filesystem::path baseName = dirPathFS.filename().empty() ? dirPathFS.parent_path().filename() : dirPathFS.filename();
	//std::cout << baseName << std::endl;
	
	// Make the output path in the present directory with the same directory name 
	std::filesystem::path outputDir = outputBaseDir / baseName;
	std::filesystem::create_directories(outputDir);
	//std::cout << outputDir << std::endl;
	
	// Spacing for the output file name
	const std::vector<double> spacing = {2*mapHalfSize[0] / nVoxels[0], 2 * mapHalfSize[1] / nVoxels[1], 2 * mapHalfSize[2] / nVoxels[2]};
	std::string mergedOutputName = gantryName + "_" + std::to_string((int)std::round(spacing[0])) + "_" + std::to_string((int)std::round(spacing[1])) + "_" + std::to_string((int)std::round(spacing[2]));
	//std::cout << mergedOutputName << std::endl;
	

	// Run sequentially
	for (unsigned int ii = 0; ii < fullPaths.size(); ii++) {
		//std::cout << fullPaths[ii] << std::endl;
		//std::cout << outputDir.string() + "/" + gSystem->BaseName(fullPaths[ii]) << std::endl;
		TH3D* h_map = getSensitivityMap(fullPaths[ii], treeName, gantryName, mapCenter, mapHalfSize, nVoxels);
		h_map->SaveAs((outputDir.string() + "/" + gSystem->BaseName(fullPaths[ii])).c_str());
		delete h_map;
	}

	
	/*
	// Run in parallel with ROOT thread pool: TTree has internal locks
	// that slow donw the I/O; this approach cannot max out the system
	auto h_maps = getSensitivityMapsThreads(fullPaths, treeName, gantryName, mapCenter, mapHalfSize, nVoxels);
	
	for (unsigned int ii = 0; ii < h_maps.size(); ii++) {
		//std::cout << fullPaths[ii] << std::endl;
		//std::cout << h_maps[ii]->GetName() << std::endl;
		h_maps[ii]->SaveAs((outputDir.string() + "/" + gSystem->BaseName(fullPaths[ii])).c_str());
		delete h_maps[ii];
	}
	h_maps.clear();
	*/


	
	// Run in parallel with separate processes
	getSensitivityMapsProcesses(fullPaths, treeName, gantryName, mapCenter, mapHalfSize, nVoxels, outputDir.string());
	
	// Merge
	mergeAllSensitivityMaps(outputDir.string(), mergedOutputName);
}


int main(int argc, char* argv[]) {
	// Parameters (todo: remove hard coding)
	/*
	const std::string simulationsBasePath = "/data/local1/raedler/J-PET/Gate_mac_9.4/New_TB_J-PET_Brain/Output/Sensitivity_brain_insert_4_18";
	const std::vector<std::string> simulationsDirs = {
		"2025-10-02_15-28-11",  // 1
		"2025-10-03_15-17-17",  // 2
		"2025-10-07_15-28-09",  // 3
		"2025-10-10_17-44-45",  // 4
		"2025-10-21_11-20-21",  // 5
		"2025-10-23_15-26-48",  // 6
		"2025-12-05_15-06-14",  // 7
		"2025-12-12_11-22-20",  // 8
		"2025-12-13_22-50-07",  // 9
		"2026-01-01_19-34-32"}; // 10
	
	// Voxelization parameters (todo: remove hard coding)
	const std::vector<double> mapCenter = {0., 0., 755.};  // mm
	const std::vector<double> mapHalfSize = {144., 144., 165.};  // mm
	//const std::vector<int> nVoxels = {288, 288, 330};  // 1×1×1 mm spacing
	//const std::vector<int> nVoxels = {144, 144, 165}; // 2×2×2 mm spacing
	const std::vector<int> nVoxels = {36, 36, 33};  // 8×8×10 mm spacing
	*/
	
	const std::string simulationsBasePath = "/data/local1/raedler/J-PET/Gate_mac_9.4/New_TB_J-PET_Brain/Output/Sensitivity_line_source";
	const std::vector<std::string> simulationsDirs = {
		"2025-12-31_15-27-48"};
	
	// Voxelization parameters (todo: remove hard coding)
	const std::vector<double> mapCenter = {0., 0., 0.};  // mm
	const std::vector<double> mapHalfSize = {0.5, 0.5, 1270.};  // mm
	const std::vector<int> nVoxels = {1, 1, 2540};  // 1×1×1 mm spacing
	
	
	// todo: remove hard coding
	const std::string treeName = "MergedCoincidences";
	
	// Choose gantry
	std::string gantryName = "Comb.";
	//std::string gantryName = "TB-TB";
	//std::string gantryName = "TB-BI";
	//std::string gantryName = "BI-BI";

	// some change
	
	// For the output directory
	std::filesystem::path simulationsBasePathFS(simulationsBasePath);
	std::filesystem::path simulationsBaseDir = simulationsBasePathFS.filename().empty() ? simulationsBasePathFS.parent_path().filename() : simulationsBasePathFS.filename();
	
	std::filesystem::path execPath = std::filesystem::canonical(argv[0]);
	std::filesystem::path execDir  = execPath.parent_path();
	
	std::filesystem::path outputBaseDir = execDir / "Output" / simulationsBaseDir;
	std::cout << outputBaseDir << std::endl;
	
	for (size_t ii = 0; ii < simulationsDirs.size(); ++ii) {
		processOneOutputDirectory(outputBaseDir, simulationsBasePath + "/" + simulationsDirs[ii] + "/", gantryName, treeName, mapCenter, mapHalfSize, nVoxels);
	}
	
	return 0;
}
