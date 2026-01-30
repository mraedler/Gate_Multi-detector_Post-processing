#ifndef CASTORTOOLS_H
#define CASTORTOOLS_H

#include <vector>
#include <utils.h>

struct LutEntry {
    float Posx;
    float Posy;
    float Posz;
    float OrVx;
    float OrVy;
    float OrVz;
};

struct ScannerParams {
    // The gaussParams and the longitudinalSpacing must be in the same length unit, e.g. [mm]
    std::vector<std::vector<std::uint32_t>> gantryShape;
    std::vector<std::vector<double>> gaussParams;
    double longitudinalSpacing;
};

std::vector<LutEntry> readLutBinary(const char* filePath);

ScannerParams totalBodyJPETWithBrainInsert_4_18();

std::vector<int> getCastorID(TTree* tree, const TLeaf* gantryID, const TLeaf* rsectorID, const TLeaf* crystalID, const TLeaf* layerID);
void checkCastorID(TTree* tree, const TLeaf* gantryID, const  TLeaf* globalPosX, const TLeaf* globalPosY, const TLeaf* globalPosZ, const std::vector<int>& castorIDs, const std::vector<LutEntry>& lut);

#endif //CASTORTOOLS_H