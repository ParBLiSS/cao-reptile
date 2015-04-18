/***
 **
 *  File: ECDriver.h
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

class ECDriver {
public:
    ECDriver(ECData& p, std::string& fn, Para& mpara):
        ecdata_(p),outFname_(fn),inPara_(mpara),readID_(0){};
    virtual ~ECDriver(){};
    void load ();
    void reLoad();
    void tableMaker();
    void ec();
private:
    ECData& ecdata_;
    std::string outFname_;
    Para& inPara_;
    int readID_;
    std::ofstream oHandle_;

    // Variables for tables
    iivec_t maskIdx_;
    uvec_t masks_;
    uvec_t table_;

    // Chunk divisonal logic taken from Reptile
    void tableUnitNeighbor(uvec_t& myNB, int dPoint);

    void processBatch(const cvec_t &ReadsString,const cvec_t &QualsString,
                      const ivec_t &ReadsOffset,const ivec_t &QualsOffset);

    void processBatchST(const cvec_t &ReadsString,const cvec_t &QualsString,
                        const ivec_t &ReadsOffset,const ivec_t &QualsOffset);
    void processBatchMT(const cvec_t &ReadsString,const cvec_t &QualsString,
                        const ivec_t &ReadsOffset,const ivec_t &QualsOffset);
    void processReadsFromFile();
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
