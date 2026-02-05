#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

#include "CASToRTools.h"
#include <TTree.h>
#include <TLeaf.h>

void identifyTrueEvents(TTree* tree, Bool_t& bBool, TBranch* b);

void setScatterTest(TTree* tree, Int_t& castorID1, Int_t& castorID2, std::vector<LutEntry>& lut, Double_t& scatterTest, TBranch* b, bool verbose);

void runMinSectorDifferenceTest(TTree* tree,
    const TLeaf* gantryID1, const TLeaf* gantryID2,
    const TLeaf* rsectorID1, const TLeaf* rsectorID2,
    Bool_t& bBool, TBranch* b,
    int minSectorDifference,
    bool verbose);

double findFirstMinimumAfterZero(TH1D* hist, bool verbose);

void runScatterTest(TTree* tree,
    const TLeaf* time1, const TLeaf* time2,
    const TLeaf* castorID1, const TLeaf* castorID2,
    const TLeaf* gantryID1, const TLeaf* gantryID2,
    const TLeaf* trueness,
    std::vector<LutEntry>& lut,
    Bool_t& bBool, TBranch* b,
    bool verbose);

void coincidenceGrouping(TTree* tree, Int_t& groupID, TBranch* b0, Int_t& groupEventID1, TBranch* b1, Int_t& groupEventID2, TBranch* b2, bool verbose = false);

std::vector<Long64_t> getGroupEdges(TTree* tree);
void selectBasedOnTime(TTree* tree, Bool_t& selection, TBranch* b, bool verbose);
void selectBasedOnEnergy(TTree* tree, Bool_t& selection, TBranch* b, bool verbose);

#endif //EVENTSELECTION_H