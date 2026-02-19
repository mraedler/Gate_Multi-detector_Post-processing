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

#pragma pack(push, 1)
struct CdfEntry {
    uint32_t time_ms;
    float    delta_time_ps;
    uint32_t castorID1;
    uint32_t castorID2;
};
#pragma pack(pop)

std::vector<LutEntry> readLutBinary(const char* filePath);

ScannerParams totalBodyJPETWithBrainInsert_4_18();
ScannerParams totalBodyJPETWithBrainInsert_6_30();

void setCastorID(TTree* tree, const TLeaf* gantryID, const TLeaf* rsectorID, const TLeaf* crystalID, const TLeaf* layerID, Int_t& castorID, TBranch* b, TString& lutName);
void checkCastorID(TTree* tree, const TLeaf* gantryID, const  TLeaf* globalPosX, const TLeaf* globalPosY, const TLeaf* globalPosZ, Int_t& castorID, const std::vector<LutEntry>& lut);

std::vector<CdfEntry> readCdfFile(const TString& filename);

#endif //CASTORTOOLS_H