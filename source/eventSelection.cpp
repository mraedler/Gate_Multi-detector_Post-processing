#include "../include/eventSelection.h"

#include <iostream>
#include <iomanip>
#include <set>
#include <algorithm>
#include <numeric>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLine.h>



void identifyTrueEvents(TTree *tree, Bool_t& bBool, TBranch* b) {

	TLeaf* eventID1 = tree->GetLeaf("eventID1");
	TLeaf* eventID2 = tree->GetLeaf("eventID2");

	TLeaf* comptonCrystal1 = tree->GetLeaf("comptonCrystal1");
	TLeaf* comptonCrystal2 = tree->GetLeaf("comptonCrystal2");
	TLeaf* rayleighCrystal1 = tree->GetLeaf("RayleighCrystal1");
	TLeaf* rayleighCrystal2 = tree->GetLeaf("RayleighCrystal2");

	TLeaf* comptonPhantom1 = tree->GetLeaf("comptonPhantom1");
	TLeaf* comptonPhantom2 = tree->GetLeaf("comptonPhantom2");
	TLeaf* rayleighPhantom1 = tree->GetLeaf("RayleighPhantom1");
	TLeaf* rayleighPhantom2 = tree->GetLeaf("RayleighPhantom2");

	Long64_t nEntries = tree->GetEntries();
	float passingPercentage = 0;
	for (Long64_t ii = 0; ii < nEntries; ii++) {
		tree->GetEntry(ii);
		bool sameEvent = eventID1->GetValue() == eventID2->GetValue();

		bool firstComptonCrystal = (comptonCrystal1->GetValue() == 1) && (comptonCrystal2->GetValue() == 1);
		bool zeroRayleighCrystal = (rayleighCrystal1->GetValue() == 0) && (rayleighCrystal2->GetValue() == 0);

		bool zeroComptonPhantom = (comptonPhantom1->GetValue() == 0) && (comptonPhantom2->GetValue() == 0);
		bool zeroRayleighPhantom = (rayleighPhantom1->GetValue() == 0) && (rayleighPhantom2->GetValue() == 0);

		bool trueCoincidence = sameEvent && firstComptonCrystal && zeroRayleighCrystal && zeroComptonPhantom && zeroRayleighPhantom;

		bBool = trueCoincidence;
		passingPercentage += bBool;
		b->Fill();
	}

	std::cout << "Percentage of true events "
		  << std::fixed << std::setprecision(2)
		  << (passingPercentage / nEntries * 100)
		  << " %.\n\n";
}



void setScatterTest(TTree* tree, Int_t& castorID1, Int_t& castorID2, std::vector<LutEntry>& lut, Double_t& scatterTest, TBranch* b, const bool verbose) {
	double speedOfLight = 2.99792458e11;  // [mm / s]

	TLeaf* time1 = tree->GetLeaf("time1");
	TLeaf* time2 = tree->GetLeaf("time2");

	TH1D* hist = nullptr;
	if (verbose) {hist = new TH1D("h", ";Scatter test [cm]; Count", 401, -100.5, 300.5);}

	Long64_t nEntries = tree->GetEntries();
	for (Long64_t ii = 0; ii < nEntries; ++ii) {
		tree->GetEntry(ii);

		int cID1 = castorID1;
		int cID2 = castorID2;

		double t1 = time1->GetValue();
		double t2 = time2->GetValue();

		float x1 = lut[cID1].Posx;
		float y1 = lut[cID1].Posy;
		float z1 = lut[cID1].Posz;

		float x2 = lut[cID2].Posx;
		float y2 = lut[cID2].Posy;
		float z2 = lut[cID2].Posz;

		double distance = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));  // [mm]
		scatterTest = (distance - (t2 - t1) * speedOfLight) / 10.;  // [cm]
		if (verbose) {hist->Fill(scatterTest);}
		b->Fill();
	}

	if (verbose) {
		double binCenterMinimum = findFirstMinimumAfterZero(hist, verbose);

		TApplication app("app", 0, nullptr);
		TCanvas canvas("c", "c", 800, 600);
		hist->Draw();
		hist->SetStats(0);

		canvas.Update();
		TLine* line = new TLine(binCenterMinimum, gPad->GetUymin(), binCenterMinimum, gPad->GetUymax());
		line->SetLineColor(kRed);
		line->SetLineWidth(2);
		line->Draw("same");

		app.Run();
	}
}



void runMinSectorDifferenceTest(TTree* tree, const TLeaf* gantryID1, const TLeaf* gantryID2, const TLeaf* rsectorID1, const TLeaf* rsectorID2, Bool_t& bBool, TBranch* b, int minSectorDifference, bool verbose) {
    // const ScannerParams sp = totalBodyJPETWithBrainInsert_4_18();
    const ScannerParams sp = totalBodyJPETWithBrainInsert_6_30();

    TH1D* hist = nullptr;
    if (verbose) {hist = new TH1D("h", ";Sector difference; Count", 15, -2.5, 12.5);}

    Long64_t nEntries = tree->GetEntries();
	double passingPercentage = 0;
    for (Long64_t ii = 0; ii < nEntries; ii++) {
        tree->GetEntry(ii);

        int gID1 = gantryID1->GetValue();
        int gID2 = gantryID2->GetValue();

        int rsID1 = rsectorID1->GetValue();
        int rsID2 = rsectorID2->GetValue();

        int sectorDifference = -1;
        if (gID1 == gID2) {
            int nSectors = sp.gantryShape[gID1][2];
            int absolute_sector_difference = std::abs(rsID2 - rsID1);
            sectorDifference = std::min(absolute_sector_difference, nSectors - absolute_sector_difference);
        }

        if (verbose) {hist->Fill(sectorDifference);}

        if (sectorDifference == -1 or sectorDifference >= minSectorDifference) {bBool = true;}
        else {bBool = false;}
		passingPercentage += bBool;
        b->Fill();
    }

	std::cout << "Passing percentage after setting a minimum sector difference of "
			  << minSectorDifference << ": "
			  << std::fixed << std::setprecision(2)
			  << (passingPercentage / nEntries * 100)
			  << " %.\n";

    if (verbose) {
        TApplication app("app", 0, nullptr);
        TCanvas canvas("c", "c", 800, 600);
        hist->Draw();
        hist->SetStats(0);
        app.Run();
    }
}



double findFirstMinimumAfterZero(TH1D* hist, bool verbose) {
	double binCenterMinimum;

	int nBins = hist->GetNbinsX();
	int jj0 = hist->FindBin(0.);

	if (jj0 == 1 || jj0 == nBins) {
		throw std::runtime_error("Zero is outside histogram range.");
	}

	for (int jj = jj0; jj <= nBins - 1; ++jj) {
		double binsCenter = (hist->GetBinCenter(jj + 1) + hist->GetBinCenter(jj)) / 2;
		double countDifference = hist->GetBinContent(jj + 1) - hist->GetBinContent(jj);
		if (verbose) {std::cout << jj << "  binsCenter = " << binsCenter << "  countDifference = " << countDifference << std::endl;}

		// Search for the point where the derivative changes
		if (countDifference > 0.) {
			binCenterMinimum = hist->GetBinCenter(jj);
			break;
		}

	}

	if (verbose) {std::cout << binCenterMinimum << std::endl;}

	return binCenterMinimum;
}



void runScatterTest(TTree* tree, const TLeaf* scatterTest, const TLeaf* gantryID1, const TLeaf* gantryID2, const TLeaf* trueness, Bool_t& bBool, TBranch* b, bool verbose) {
	Long64_t nEntries = tree->GetEntries();

	// Different thresholds for different gantries
	// todo: remove hard coding
	double thresholdTBTB = 25.;  // [cm]
	double thresholdOther = 10.;  // [cm]

	TH1D* hist = nullptr;
	if (verbose) {hist = new TH1D("h", ";Scatter test [cm]; Count", 401, -100.5, 300.5);}

	double passingPercentage = 0;
	double excludedTrue = 0;
	double trueCounter = 0;
	for (Long64_t ii = 0; ii < nEntries; ++ii) {
		tree->GetEntry(ii);

		double st = scatterTest->GetValue();

		// Only include the entries that were not already excluded by the minimum sector difference
		if (verbose) {if (bBool) hist->Fill(st);}

		if ((gantryID1->GetValue() < 2) && (gantryID2->GetValue() < 2)) {
			if (st < thresholdTBTB) {bBool = false;}
		} else {
			if (st < thresholdOther) {bBool = false;}
		}

		passingPercentage += bBool;

		bool t = trueness->GetValue();
		trueCounter += t;

		if (!bBool && t) {excludedTrue++;}
	}

	std::cout << "Passing percentage after additional scatter test: "
		  << std::fixed << std::setprecision(2)
		  << (passingPercentage / nEntries * 100)
		  << " %.\n";

	std::cout << "Percentage of excluded true events: "
	  << std::fixed << std::setprecision(2)
	  << (excludedTrue / trueCounter * 100)
	  << " %.\n\n";


	if (verbose) {
		double binCenterMinimum = findFirstMinimumAfterZero(hist, verbose);

		TApplication app("app", 0, nullptr);
		TCanvas canvas("c", "c", 800, 600);
		hist->Draw();
		hist->SetStats(0);

		canvas.Update();
		TLine* line = new TLine(binCenterMinimum, gPad->GetUymin(), binCenterMinimum, gPad->GetUymax());
		line->SetLineColor(kRed);
		line->SetLineWidth(2);
		line->Draw("same");

		app.Run();
	}


	// std::vector<bool> passingScatterTest(nEntries, false);
	// float passingCounter = 0.;
	//
	// for (Long64_t ii = 0; ii < nEntries; ++ii) {
	// 	tree->GetEntry(ii);
	//
	// 	if ((gantryID1->GetValue() < 2) && (gantryID2->GetValue() < 2)) {
	// 		if (scatterTest[ii] >= thresholdTBTB) {
	// 			passingScatterTest[ii] = true;
	// 			passingCounter++;
	// 		}
	// 	} else {
	// 		if (scatterTest[ii] >= thresholdOther) {
	// 			passingScatterTest[ii] = true;
	// 			passingCounter++;
	// 		}
	// 	}
	// }

	//if (verbose) {std::cout << "Passing percentage scatter test: " << passingCounter / nEntries * 100 << " %." << std::endl;}
}



void coincidenceGrouping(TTree* tree, Int_t& groupID, TBranch* b0, Int_t& groupEventID1, TBranch* b1, Int_t& groupEventID2, TBranch* b2, const bool verbose)
{
	TLeaf* time1 = tree->GetLeaf("time1");
	TLeaf* time2 = tree->GetLeaf("time2");
	TLeaf* energy1 = tree->GetLeaf("energy1");
	TLeaf* energy2 = tree->GetLeaf("energy2");

	// std::vector<Long64_t> groupEdges = {0};
	// std::vector<Long64_t> groupMultiplicities;

	TH1D* histMultiplicity = nullptr;
	if (verbose) {histMultiplicity = new TH1D("h1", ";Multiplicity; Count", 12, -0.5, 11.5);}

	TH1D* histGroupSize = nullptr;
	if (verbose) {histGroupSize = new TH1D("h2", ";Group size; Count", 21, -0.5, 20.5);}

	Long64_t ii = 0;
	Long64_t nEntries = tree->GetEntries();
	Int_t groupIDTemp = 0;

	while (ii < nEntries) {

		//std::set<std::tuple<double, double>> seen;
		std::map<std::tuple<double, float>, int> seen;
		int nextId = 0;
		Int_t groupSize = 0;
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

			groupEventID1 = seen[k1];
			groupEventID2 = seen[k2];
			b1->Fill();
			b2->Fill();

			groupID = groupIDTemp;
			b0->Fill();

			groupSize++;
			ii++;
		}

		groupIDTemp++;
		if (verbose) {histMultiplicity->Fill(seen.size());}
		if (verbose) {histGroupSize->Fill(groupSize);}
		// groupEdges.push_back(ii);
		// groupMultiplicities.push_back(seen.size());
	}

	/*
	// Reconstruct the groupEdges from the groupID
	std::vector<Long64_t> groupEdges = {0};
	Long64_t previousGroupID = 0;
	for (Long64_t jj = 1; jj < nEntries; ++jj) {
		tree->GetEntry(jj);
		if (groupID - previousGroupID > 0) {
			groupEdges.push_back(jj);
		}
		previousGroupID = groupID;
	}
	groupEdges.push_back(nEntries);

	if (groupEdges.back() != nEntries) {
	std::cerr << "Warning: not all elements grouped!\n";
	}
	*/

	/*
	// Print all time entries
	for (Long64_t jj = 0; jj < nEntries; ++jj) {
		tree->GetEntry(jj);
		std::cout << std::setprecision(16) << jj << " " << time1->GetValue() << " " << time2->GetValue() << std::endl;
	}
	std::cout << std::endl;
	*/

	/*
	// Print grouped time entries
	for (Long64_t jj = groupEdges.size() - 100; jj < groupEdges.size() - 1; ++jj){
		std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;
		for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; ++kk) {
			tree->GetEntry(kk);
			std::cout << std::setprecision(16) << time1->GetValue() << " " << time2->GetValue() << std::endl;
		}
		std::cout << std::endl;
	}
	*/

	/*
	// Check for multiplicity 1, indicating whether the criterion (time only) to identify identical events is insufficient
	for (Long64_t jj = 0; jj < groupEdges.size() - 1; ++jj){
		if (groupMultiplicities[jj] == 1) {
			std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;
			for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; ++kk) {
				tree->GetEntry(kk);
				//std::cout << std::setprecision(16) << time1->GetValue() << " " << time2->GetValue() << std::endl;
				std::cout << std::setprecision(16) << time1->GetValue() - time2->GetValue() << std::endl;
				std::cout << std::setprecision(16) << energy1->GetValue() - energy2->GetValue() << std::endl;
			}
			std::cout << std::endl;
		}
	}
	*/

	/*
	// Check the indexing
	for (Long64_t jj = 0; jj < groupEdges.size() - 1; ++jj){
		std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;
		for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; ++kk) {
			tree->GetEntry(kk);
			std::cout << groupID << " " << std::setprecision(16) << time1->GetValue() << " " << groupEventID1 << " " << time2->GetValue() << " " << groupEventID2 << std::endl;
		}
		std::cout << std::endl;
	}
	*/

	if (verbose) {
		TApplication app("app", 0, nullptr);
		TCanvas canvas("c", "c", 1200, 600);
		canvas.Divide(2,1);

		// Expected to peak at the triangle numbers for the data before the pre-processing: 1, 3, 6, 10, 15, 21, ...
		canvas.cd(1);
		gPad->SetLogy();
		histGroupSize->Draw();
		histGroupSize->SetStats(0);

		canvas.cd(2);
		gPad->SetLogy();
		histMultiplicity->Draw();
		histMultiplicity->SetStats(0);

		app.Run();
	}
}



std::vector<Long64_t> getGroupEdges(TTree* tree) {
	Long64_t nEntries = tree->GetEntries();
	TLeaf* groupID = tree->GetLeaf("groupID");

	std::vector<Long64_t> groupEdges;
	groupEdges.reserve(nEntries);
	groupEdges.push_back(0);

	Long64_t previousGroupID = 0;
	for (Long64_t jj = 1; jj < nEntries; ++jj) {
		tree->GetEntry(jj);
		if (groupID->GetValue() - previousGroupID > 0) {
			groupEdges.push_back(jj);
		}
		previousGroupID = groupID->GetValue();
	}

	groupEdges.push_back(nEntries);

	return groupEdges;
}



void selectBasedOnTime(TTree* tree, Bool_t& selection, TBranch* b, const bool verbose) {

	const std::vector<Long64_t> groupEdges = getGroupEdges(tree);

	TLeaf* groupEventID1 = tree->GetLeaf("groupEventID1");
	TLeaf* groupEventID2 = tree->GetLeaf("groupEventID2");

	// std::vector<bool> selection(idx1.size(), false);

	for (Long64_t jj = 0; jj < groupEdges.size() - 1; jj++){
		if (verbose) {std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;}

		std::set<Long64_t> usedEvents;
		for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; kk++) {
			tree->GetEntry(kk);
			if (verbose) {std::cout << groupEventID1->GetValue() << " " << groupEventID2->GetValue() << " ";}

			// if (usedEvents.count(groupEventID1->GetValue()) || usedEvents.count(groupEventID2->GetValue())) {
			if (usedEvents.count(groupEventID1->GetValue()) || usedEvents.count(groupEventID2->GetValue()) || usedEvents.size() > 0) {  // only choose one coincidence
				if (verbose) {std::cout << std::endl;}
				selection = false;
				b->Fill();
				continue;
			}

			if (verbose) {std::cout << "*" << std::endl;}

			selection = true;
			b->Fill();
			usedEvents.insert(groupEventID1->GetValue());
			usedEvents.insert(groupEventID2->GetValue());

		}
		if (verbose) {std::cout << std::endl;}
	}
}



void selectBasedOnEnergy(TTree* tree, Bool_t& selection, TBranch* b, const bool verbose) {

	const std::vector<Long64_t> groupEdges = getGroupEdges(tree);

	TLeaf* energy1 = tree->GetLeaf("energy1");
	TLeaf* energy2 = tree->GetLeaf("energy2");

	TLeaf* groupEventID1 = tree->GetLeaf("groupEventID1");
	TLeaf* groupEventID2 = tree->GetLeaf("groupEventID2");

	for (Long64_t jj = 0; jj < groupEdges.size() - 1; jj++){
		if (verbose) {std::cout << groupEdges[jj] << " " << groupEdges[jj + 1] << std::endl;}

		// Cache total energy
		std::vector<float> totalEnergies;
		totalEnergies.reserve(groupEdges[jj+1] - groupEdges[jj]);
		for (Long64_t kk = groupEdges[jj]; kk < groupEdges[jj + 1]; ++kk) {
			tree->GetEntry(kk);
			totalEnergies.push_back(energy1->GetValue() + energy2->GetValue());
		}

		// Arg-sort based on the energy
		std::vector<Long64_t> indices(groupEdges[jj + 1] - groupEdges[jj]);
		std::iota(indices.begin(), indices.end(), 0);
		std::sort(indices.begin(), indices.end(), [&](std::size_t aa, std::size_t bb) {return totalEnergies[aa] > totalEnergies[bb];});

		// The branch can only be written sequentially; first need to figure out the order
		std::vector<Bool_t> selectionGroup(groupEdges[jj + 1] - groupEdges[jj]);

		// Processing out of order
		std::set<Long64_t> usedEvents;
		for (Long64_t kk = 0; kk < indices.size(); kk++) {
			// int sortedGroupIndex = sortedGroupIndices[ll];
			tree->GetEntry(groupEdges[jj] + indices[kk]);
			if (verbose) {std::cout << groupEventID1->GetValue() << " " << groupEventID2->GetValue() << " " << energy1->GetValue() + energy2->GetValue() << " ";}

			//if (usedEvents.count(groupEventID1->GetValue()) || usedEvents.count(groupEventID2->GetValue())) {
			if (usedEvents.count(groupEventID1->GetValue()) || usedEvents.count(groupEventID2->GetValue()) || usedEvents.size() > 0) {  // only choose one coincidence
				if (verbose) {std::cout << std::endl;}
				selectionGroup[indices[kk]] = false;
				// selection = false;
				// b->Fill();
				continue;
			}

			if (verbose) {std::cout << "*" << std::endl;}

			selectionGroup[indices[kk]] = true;
			// selection = true;
			// b->Fill();
			usedEvents.insert(groupEventID1->GetValue());
			usedEvents.insert(groupEventID2->GetValue());
		}
		if (verbose) {std::cout << std::endl;}

		// Fill in order
		for (Long64_t kk = 0; kk < selectionGroup.size(); kk++) {
			selection = selectionGroup[kk];
			b->Fill();
		}
	}
}