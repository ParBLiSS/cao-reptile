/***
 **
 *  File: Parser.h
 *  Created: Dec 12, 2009
 *
 *  Author: Xiao Yang <isuxyang@gmail.com>
 *
 *  Copyright 2009 Xiao Yang, Karin Dorman, Srinivas Aluru
 *  Copyright 2011-12 Ankit Shah, Sriram P C, Srinivas Aluru
 *
 *  This file is part of Reptile (version 1.1).
 *
 *  Reptile is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Reptile is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _PARSER_H
#define	_PARSER_H

#include "util.h"
#include "fasta_file.hpp"
#include "ECData.hpp"


typedef std::vector<kc_t> kcvec_t;

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


class Parser {
public:
    Parser(ECData *p){ readID_ = 0; ecdata = p;};
    Parser(const Parser& orig){};
    virtual ~Parser(){};
    void load (const Para& myPara);
    void reLoad(const Para& myPara);
    void tableMaker(const Para& myPara);
    void ec(const Para& myPara);
    void output (const std::string& filename);
private:
    std::vector<record_t> records_;
    ECData *ecdata;

    // Variables for tables
    iivec_t maskIdx_;
    uvec_t masks_;
    uvec_t table_;

    // temporary variables
    int readID_;
    evec_t readErr_;    

    void readEC(char* addr, char* qAddr, const Para& myPara);
    bool errorCall(const upair_t& mosaic, const ipair_t& dPoints,
                const ipair_t& hdUb, char* qAddr, const Para& myPara);

    // Chunk divisonal logic taken from Reptile
    void tableMaker();
    void tableUnitNeighbor(uvec_t& myNB, int dPoint,const Para& myPara);

    // Find neighbor searching in all kmer list
    void genericUnitNeighbor(uvec_t& myNB, int dPoint,const Para& myPara);

    void unitNeighbor(uvec_t& myNB, int dPoint, const Para& myPara);
    void updateRead(std::string& myRead, int shift);
    void neighbors(std::vector<kcvec_t>& nodes,
        const uvec_t& repIDs, const Para& myPara);
    void errCall (std::vector<kcvec_t>& nodes, std::vector<kcvec_t>& nodes_rv,
        const std::string& qual, const Para& myPara);
    void tiling(kcvec_t& tiles, const uvec_t& N1, const uvec_t& N2, const Para& myPara);
    void mergeTiles(kcvec_t& tileTo, kcvec_t& tileFrom, const Para& myPara);
    bool overlay(uint64_t& rslt, uint32_t n1, uint32_t n2, const Para& myPara);

    void processBatch(cvec_t &ReadsString,cvec_t &QualsString,
                      ivec_t &ReadsOffset,ivec_t &QualsOffset,
                      const Para& myPara);
    void processReadsFromFile(const Para& myPara);
};

void updateNodes(kcvec_t& N, kcvec_t& N_rv, const kcvec_t& tiles, int len);

// used for sorting tables
struct TComp {
    uint32_t mask;
    TComp (uint32_t mvalue): mask(mvalue) {};

    bool operator() (uint32_t e1, uint32_t e2) const {
        return ((e1 & mask) < (e2 & mask));
    }
};

/*
 * for debugging purpose, print out std::vector<kc_t>
 */
inline void print_kcvec (const std::string& msg, const kcvec_t& myvec, int len){

    std::cout << msg << "\n";
    std::cout << "\n--------------------------------\n";
    for (unsigned int j = 0; j < myvec.size(); ++ j){
        std::cout << toString(myvec[j].ID, len)
        << "\t" << myvec[j].goodCnt << "\t" << myvec[j].cnt << "\n";
    }
    std::cout << "--------------------------------\n";
}

#endif	/* _PARSER_H */


