#ifndef ERRORCORRECTOR_H
#define ERRORCORRECTOR_H

#include <vector>
#include <string>
#include <ostream>
#include "util.h"
#include "ECData.hpp"

typedef struct ERRINFO{
    int pos;
    int from;
    int to;
    int qual;
    ERRINFO (int p, int f, int t, int q): pos(p), from(f), to(t), qual(q) {}
    ERRINFO() {};
}e_t;
typedef std::vector<e_t> evec_t;

typedef struct RECORD{
    int readID;
    evec_t evec;
    RECORD (int ID, evec_t vec): readID(ID), evec(vec) {}
}record_t;


class ECImpl{
private:
    std::vector<record_t> records_;
    ECData& ecdata_;
    Para& inPara_;
    // temporary variables
    evec_t readErr_;

    bool errorCall(const upair_t& mosaic, const ipair_t& dPoints,
                   const ipair_t& hdUb, char* qAddr);
    void updateRead(int readID, std::string& myRead, int shift);
    bool correctTile(const kcvec_t& tiles, const char* qAddr);
    void tiling(kcvec_t& tiles, const uvec_t& N1, const uvec_t& N2);
    bool overlay(uint64_t& rslt, uint32_t n1, uint32_t n2);

    void mergeTiles(kcvec_t& tileTo, kcvec_t& tileFrom);
    void genericUnitNeighbor(uvec_t& myNB, int dPoint);
    void unitNeighbor(uvec_t& myNB, int dPoint);

public:
    ECImpl(ECData& ecd, Para& inp):ecdata_(ecd), inPara_(inp){};
    void readEC(int readID, char* addr, char* qAddr);
    void writeErrors(std::ostream& oHandle);
};

#endif /* ERRORCORRECTOR_H */
