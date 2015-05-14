#include <cmath>

#include "util.h"
#include "find_neighbors.h"
#include "ECImpl.hpp"

void candidates(ivec_t& candies, const kcvec_t& tiles, int threshold) {
    for (unsigned i = 1; i < tiles.size(); ++i) {
        if (tiles[i].goodCnt >= threshold)
            candies.push_back(i);
    }
}

void diff(evec_t& errs, uint64_t from, uint64_t to, const char* qAddr, int len) {
    uint64_t tmp = from ^ to;
    int idx = 0;
    while (tmp) {
        if (tmp & 0x3) {
            e_t err;
            err.pos = len - 1 - idx;
            err.from = from & 0x3;
            err.to = to & 0x3;
            err.qual = (int) qAddr[(len - 1) - idx];
            errs.push_back(err);
        }
        tmp >>= 2;
        from >>= 2;
        to >>= 2;
        idx++;
    }
}

bool goodQuals(const char* qAddr, int len, int threshold) {
    for (int i = 0; i < len; ++i) {
        if (qAddr[i] < threshold)
            return false;
    }
    return true;
}

bool unicpy(const kc_t& e1, const kc_t& e2) {
    return (e1.ID == e2.ID);
}

void ECImpl::readEC(int readID, char* addr, char* qAddr) {


    ipair_t hdUb(inPara_.hdMax, inPara_.hdMax);

    ipair_t dPoints(0, inPara_.K - inPara_.step); //dividing points //a

    std::string myRead = addr;
    int readLen = myRead.length();
    int nextPos = 0, prePos = -1;
    bool crawlingFlag = false;
    int crawlTo = -1;

    while (1) {

        upair_t mosaic;
        uint64_t ID1, ID2;
        if (!(toID(ID1, &myRead[nextPos], inPara_.K) &&
              toID(ID2, &myRead[nextPos + inPara_.step], inPara_.K))) {
            return; /* containing non-acgt*/
        }
        mosaic = upair_t(ID1, ID2);

        if (errorCall(mosaic, dPoints, hdUb, qAddr + nextPos)) {//a

            updateRead(readID, myRead, nextPos);

            if (nextPos + inPara_.step + inPara_.K == readLen) {
                return;
            }

            prePos = nextPos;

            if (!crawlingFlag) {
                if (nextPos + inPara_.step <= readLen - inPara_.K - inPara_.step) {
                    dPoints = ipair_t(inPara_.K, 0);
                    nextPos += inPara_.step;
                    hdUb = ipair_t(0, inPara_.hdMax);
                }
                else {
                    dPoints = ipair_t(inPara_.K, nextPos + inPara_.K + 2 * inPara_.step - readLen);
                    nextPos = readLen - inPara_.K - inPara_.step;
                    hdUb = ipair_t(0, inPara_.hdMax);
                }
            }
        }
        else {
            if (crawlingFlag == false){
                crawlingFlag = true;
                crawlTo = nextPos;
            }
        }

        if (crawlingFlag){

         // crawling --- to do reverse crawling
            if (nextPos == prePos + 1 || prePos == -1){
                return;
            }

            nextPos = prePos + 1;

            if(nextPos == crawlTo) crawlingFlag = false;
            hdUb = ipair_t(0, 1);
            dPoints = ipair_t(inPara_.K, inPara_.K - 3);
        }
    }
}

bool ECImpl::errorCall(const upair_t& mosaic, const ipair_t& dPoints,
        const ipair_t& hdUb, char* qAddr) {

    upair_t rvMosaic(
            reverse_complementary<uint32_t, uint32_t > (mosaic.first, inPara_.K),
            reverse_complementary<uint32_t, uint32_t > (mosaic.second, inPara_.K));

    uvec_t heads, tails, rvHeads, rvTails;
    if (ecdata_.findKmer((kmer_id_t)mosaic.first)){
        heads.push_back(mosaic.first);
    }

    if (ecdata_.findKmer((kmer_id_t)mosaic.second)){
        tails.push_back(mosaic.second);
    }

    if (ecdata_.findKmer((kmer_id_t)rvMosaic.second)){
        rvHeads.push_back(rvMosaic.second);
    }

    if (ecdata_.findKmer((kmer_id_t)rvMosaic.first)){
        rvTails.push_back(rvMosaic.first);
    }

    ipair_t hd(0, 0);
    while (hd.first <= hdUb.first) {

        kcvec_t tiles, rvTiles;
        tiling(tiles, heads, tails);
  //      print_kcvec("tiles", tiles, inPara_.K + inPara_.step);
        tiling(rvTiles, rvHeads, rvTails);
  //      print_kcvec("rvtiles", rvTiles, inPara_.K + inPara_.step);
        mergeTiles(tiles, rvTiles);
  //      print_kcvec("After Merging", tiles, inPara_.K + inPara_.step);

        if(correctTile(tiles, qAddr))
            return true;

        /*
         * increase Hamming Distance to search for more neighbors
         */
        hd.second++;
        if (hd.second <= hdUb.second) {
            unitNeighbor(tails, dPoints.second);
            unitNeighbor(rvHeads, dPoints.second);
        } else {
            hd.first++;
            if (hd.first <= hdUb.first) {
                unitNeighbor(heads, dPoints.first);
                unitNeighbor(rvTails, dPoints.first);
            }
            else{ // all possible searches have been done, till this stage
                  // if tile[0] is the maximum and > tCard then, it is considered
                  // as correct due to low coverage region
                int maxCnt = 0;
                bool tflag = true;
                if (tiles.size() > 0) maxCnt = tiles[0].goodCnt;
                if (maxCnt >= inPara_.tCard) {
                    for (unsigned i = 1; i < tiles.size(); i ++){
                        if (maxCnt < tiles[i].goodCnt) {
                            tflag = false;
                            break;
                        }
                    }
                    if (tflag) return true;
                }
            }
        }
    }
    return false;
}

void ECImpl::updateRead(int readID, std::string& myRead, int shift) {

    if (readErr_.size()) {
        // 1. adding shift to err positions w.r.t read & update read
        for (unsigned int i = 0; i < readErr_.size(); ++i) {
            readErr_[i].pos += shift;
            myRead[readErr_[i].pos] = bits_to_char(readErr_[i].from);
        }
        // 2. update records_
        if (records_.size() && records_.back().readID == readID) {
            (records_.back()).evec.insert((records_.back()).evec.end(),
                    readErr_.begin(), readErr_.end());
        } else {
            records_.push_back(record_t(readID, readErr_));
        }

        // 3. update snKArray_ and lnKArray_ --- todo

        readErr_.clear();
    }
}

bool ECImpl::correctTile(const kcvec_t& tiles, const char* qAddr) {
    /*
     * Error Calling Current Tile
     * tiles.size() may equal to 0 due to error correction
     */
    if (tiles.size() == 0){
        return false;
    }

    if (tiles[0].goodCnt >= inPara_.tGoodTile) {
        return true;
    }

    if (tiles.size() == 1) { // try to correct
        if (tiles[0].goodCnt >= inPara_.tCard &&
            goodQuals(qAddr, inPara_.K + inPara_.step, inPara_.Qlb)) {
            return true;
        }
    } else {
        ivec_t candies; // indices of tiles
        candidates(candies, tiles, inPara_.tCard);

        if (tiles[0].goodCnt >= inPara_.tCard) {
            /*
             * ec only if non-ambig correction of low qual bases
             * could be identified and correcting to one of the
             * high cardinality neighbors
             */
            int tGoodCnt = tiles[0].goodCnt / inPara_.tRatio;

            ivec_t highCardNbs;
            for (unsigned i = 0; i < candies.size(); ++i) {
                if (tiles[candies[i]].goodCnt >= tGoodCnt)
                    highCardNbs.push_back(candies[i]);
            }
            int alterNum = 0;
            for (unsigned i = 0; i < highCardNbs.size(); ++i) {
                /*
                 * check all candidates, if ambig, then do not ec
                 */
                evec_t errs;
                diff(errs, tiles[highCardNbs[i]].ID, tiles[0].ID,
                     qAddr, inPara_.step + inPara_.K);
                for (unsigned j = 0; j < errs.size(); ++j) {
                    if (errs[j].qual < inPara_.Qlb) {
                        readErr_ = errs;
                        alterNum++;
                        break;
                    }
                }
            }
            if (alterNum > 1) { // ambig
                readErr_.clear();
            } else
                return true;
        } else {
            /*
             * ec only if non-ambig correction to one of [candies]
             * could be identified
             */
            if (candies.size() == 1) {
                diff(readErr_, tiles[candies[0]].ID, tiles[0].ID,
                     qAddr, inPara_.step + inPara_.K);
                return true;
            }
        }
    }

    return false;
}

/*
 * flag: true -- forward    flase--reverse
 */
void ECImpl::tiling(kcvec_t& tiles, const uvec_t& N1,
        const uvec_t& N2) {

    if (N1.size() == 0 || N2.size() == 0) return;
    uint64_t reptile;
    if (!overlay(reptile, N1[0], N2[0])) {
        std::cout << "Err: errCall, reptile construction fail\n";
        exit(1);
    }
    kc_t output;
    long idx = ecdata_.findTile(reptile,output);
    if (idx != -1) {
        tiles.push_back(output);
    }
    else {
        tiles.push_back(kc_t(reptile, 0, 0));
    }
    for (unsigned i = 0; i < N1.size(); ++i) {
        for (unsigned j = 0; j < N2.size(); ++j) {
            if (i == 0 && j == 0) continue;
            uint64_t tmptile;
            if (overlay(tmptile, N1[i], N2[j])) {
                //binary search
                kc_t tmpoutput;
                idx = ecdata_.findTile(tmptile,tmpoutput);

                if (idx != -1) {
                    tiles.push_back(tmpoutput);
                }
            }
        }
    }
}

bool ECImpl::overlay(uint64_t& rslt, uint32_t n1, uint32_t n2) {

    int shift = 2 * inPara_.step;
    uint64_t a = n1, b = n2;
    // comparison suffix - prefix if necessary
    if (inPara_.step >= inPara_.K ||
            (b >> shift) == (a & ((1 << 2*(inPara_.K - inPara_.step)) - 1))) {
        rslt = (a << shift) | b;
        return true;
    }
    return false;
}

void ECImpl::mergeTiles(kcvec_t& tileTo, kcvec_t& tileFrom) {
    /*
     *  convert IDs in [tileFrom] to revcompl IDs and insert to [tileTo]
     *  then sort and merge (remove duplicates); keep the repID in place
     */
    if (tileFrom.size() == 0) return;

    /* convert to reverse compl IDs
     */
    for (unsigned i = 0; i < tileFrom.size(); ++i) {
        tileFrom[i].ID = reverse_complementary <uint64_t, uint64_t >
                (tileFrom[i].ID, inPara_.K + inPara_.step);

    }

    if (tileTo.size() == 0){
        tileTo = tileFrom;
        tileFrom.clear();
        return;
    }

    uint64_t repID = tileTo[0].ID;
    tileTo.insert(tileTo.end(), tileFrom.begin(), tileFrom.end());
    tileFrom.clear();

    /* sort
     */
    std::sort(tileTo.begin(), tileTo.end(), Knumcomp());

    /* merge
     */
    int idx1 = 0;

    for (int idx2 = idx1 + 1; idx2 < (int)tileTo.size(); ++idx2) {
        if (tileTo[idx2].ID == tileTo[idx1].ID) {
            tileTo[idx1].goodCnt += tileTo[idx2].goodCnt;
            tileTo[idx1].cnt += tileTo[idx2].cnt;
        } else {
            idx1 = idx2;
        }
    }

    kcvec_t::iterator it =
            std::unique_copy(tileTo.begin(), tileTo.end(), tileTo.begin(), unicpy);
    tileTo.resize(it - tileTo.begin());

    /* keep repID in the first pos
     */
    if (tileTo[0].ID != repID) {
        for (unsigned i = 1; i < tileTo.size(); ++i) {
            if (tileTo[i].ID == repID) {
                kc_t tmp = tileTo[0];
                tileTo[0] = tileTo[i];
                tileTo[i] = tmp;
                break;
            }
        }
    }
}

void ECImpl::unitNeighbor(uvec_t& myNB, int dPoint) {

    // Commented as it is no longer used
    // if(inPara_.useMaskedLists){
    //     tableUnitNeighbor(myNB, dPoint,inPara_);
    //     return;
    // }
    genericUnitNeighbor(myNB,dPoint);
}

//
// Finds the unit neighbor by generating each possible neighbor
void ECImpl::genericUnitNeighbor(uvec_t& myNB, int dPoint){

    // Parallel Reptile ::
    //    Replaced this function with the query mechanism

#ifdef DEBUG
    std::cout << "DPOINT : " <<  dPoint
              << " INPUT : " << myNB.size() << " : " ;
    for(int i = 0; i < myNB.size(); i++ )
        std::cout << myNB[i] << " " ;
    std::cout << std::endl;
#endif
    if (myNB.size() == 0 || dPoint >= inPara_.K) {
#ifdef DEBUG
        std::cout << "OUTPUT : NONE" << std::endl ;
#endif
        return;
    }

    std::vector<kmer_id_t> nbIDs;

    int unitSize = (log2(inPara_.eSearch) / 2);
    int stopPoint = -1 ;
    stopPoint =  ((dPoint+1)/unitSize)*unitSize -1;
    // TODO: update as required
    stopPoint = dPoint;
    for (unsigned i = 0; i < myNB.size(); ++i) {
        get_dneighbors(inPara_.K,1, (kmer_id_t)myNB[i],stopPoint,nbIDs);
    }

    std::sort(nbIDs.begin(), nbIDs.end());
    uvec_t::iterator it
            = std::unique_copy(nbIDs.begin(), nbIDs.end(), nbIDs.begin());
    nbIDs.resize(it - nbIDs.begin());

#ifdef DEBUG
    std::cout << "OUTPUT : " ;
#endif
    for (unsigned i = 0; i < nbIDs.size(); ++i) {
        if ( nbIDs[i] != 0 && ecdata_.findKmerCacheAware(nbIDs[i])){
#ifdef DEBUG
            std::cout << nbIDs[i] << " " ;
#endif
            myNB.push_back(nbIDs[i]);
        }
    }
#ifdef DEBUG
    std::cout << std::endl;
#endif
}

void ECImpl::writeErrors(std::ostream& oHandle) {
    /*
     * Resulting file format:
     * ReadID (sorted)  ErrNum  [pos from to qual] [pos from to qual]...
     * from: reference; to: read (numercial value Aa:0 Cc:1 Gg:2 Tt:3 others 4)
     * qual: quality (numerical value)
     */
    for (unsigned i = 0; i < records_.size(); ++i) {

        oHandle << records_[i].readID << "\t" << records_[i].evec.size();

        for (unsigned j = 0; j < records_[i].evec.size(); ++j) {
            oHandle << "\t" << records_[i].evec[j].pos << "\t"
                    << records_[i].evec[j].from << "\t" << records_[i].evec[j].to
                    << "\t" << records_[i].evec[j].qual;
        }
        oHandle << "\n";
    }

}
