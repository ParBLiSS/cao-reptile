/***
 *
 *  Author: Ankit Shah <shah.29ankit@gmail.com>
 *          Sriram P C <sriram.pc@gmail.com>
 *
 *  Copyright 2011-2012 Ankit Shah, Sriram P C, Srinivas Aluru
 *
 *  This file is part of Parallel Reptile (version 1.1)
 *
 *  PReptile is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PReptile is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "util.h"
#include <iostream>
#include <stdint.h>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>



//
// Iterating method for next r combination 
//  in the reverse order
bool next_rcombination_rev(int n, int r, int *&next) {
    int i = r-1;

    // find the next gap to decrement
    while (i > 0 && (next[i] - next[i-1] == 1) ) 
        i = i - 1;

    // reached the last one
    if(i == 0 && next[i] == 0)
        return false;

    next[i] = next[i] - 1;
    for(int j = r-1,s = n-1; j > i;j--,s--)
        next[j] = s;

    return true;
}

//
// Iterating method for next r combination
//  in the forward order
bool next_rcombination(int n,int r,int *&next){
    int i = r-1;

    // me = Find the next position to increment
    while (i >= 0 && next[i] == n - r + i ) 
        i = i - 1;

    if(i < 0)
        return false;

    // increment me and every one to my right
    next[i] = next[i] + 1;
    for(int j = i+1; j < r; j++)
        next[j] = next[i] + j - i;

    return true;
}

bool next_rcombination(int n,int r,int **pnext,int stopPoint){

    assert(r == 1); // Right now this works only if r == 1
    int i = r-1;
    int *next = *pnext;

    if(next[i] == stopPoint)
        return false;
    next[i] = next[i] - 1;
    return true;
}

template <typename T>
inline T get_neighbors(T& ID,int len,int d,int *dpositions,T nbrdata){
    static T allones = 0;
    int i = 0;

    if(!allones)
        for( i = 0; i < len; ++i)
            allones = (allones << 2) | 3;


    // nbrclear sets '11' for all the places where
    // we have to set the neighbor bits
    T nbrclear = 0;

    // nbrmask has the nbr bits on the places given by
    //  dpositions and all ones in other places
    T nbrmask = allones;

    // for the given d-positions
    //  - set the nbrdata
    for(i= 0; i < d;i++) {
        int kpos = dpositions[i];
        // nbrclear sets '11' for all the places where
        // there is kpos
        nbrclear |= (3 << (2*kpos));

        // set all ones to the right of kpos
        T nbrm = allones << (2 * (kpos+1));
        nbrm &= allones;

        // get the two bits at  position 'i' to the 
        // end of the bit vector
        T nbrbits = nbrdata >> (2*i);
        nbrbits &= 3;

        // Shift (2*kpos) and set all ones on the left of kpos
        nbrbits <<=  (2*kpos);
        nbrbits |= (allones >> (2*(len-kpos)));

        nbrmask &= (nbrbits|nbrm);
    }

#ifdef DEBUG
    std::cout << nbrclear << " " << nbrmask << " " 
              << ((ID | nbrclear) & nbrmask) << std::endl;
#endif
    return (ID | nbrclear) & nbrmask;
}

template <typename T>
void print_bits(T &ID,int l){
    std::string bitstring;
    T temp = ID;
    for(int i = 0;i < l;i++) {
        bitstring += (temp & 1) == 1 ? "1" : "0";
        temp >>= 1;
    }
    std::reverse(bitstring.begin(),bitstring.end());

    std::cout << ID << ":" << bitstring << std::endl;
}

void print_comb(int l,int *next){
    std::cout << "{" ;
    for(int j = 0; j < l;j++) 
        std::cout << next[j] << (j < l-1 ? "," :"");
    std::cout << "}" << std::endl;
}

template <typename T>
void get_neighbors(std::vector<T> &nbrs,T &ID,int len,
                   int d,int *dpositions,T &nbrdata){
    kmer_id_t max = 1 << (2*d),j = 0;
    // enumerate all possible 4^d 
    for(j = 0; j < max;j++) {
        kmer_id_t nbr = get_neighbors(ID,len,d,dpositions,j);
        if(nbr != ID) {
            nbrs.push_back(nbr);
        }
    }
}

int get_dneighbors(int k,int d,kmer_id_t ID,
         int stopPoint,std::vector<kmer_id_t> &neighbors){
    int *a = (int*) calloc(k,sizeof(int));
    int *next = (int*) a;
    int **pnext = &next;

    // for each possible value for d
    for(int i = 1; i <= d;i++) {
        kmer_id_t j = 0;
        // chose a combination from d
        // Start from the last
        // 1st d combination is just first d of 'a'
        for(int r = 0; r < d;r++) a[r] = k-r-1;
#ifdef DEBUG
        print_comb(i,next);
#endif
        get_neighbors<kmer_id_t>(neighbors,ID,k,i,next,j);
        // next d starts here
        while(next_rcombination(k,i,pnext,stopPoint)){
#ifdef DEBUG
            print_comb(i,next);
#endif
            get_neighbors(neighbors,ID,k,i,next,j);
        }
    }

    std::sort(neighbors.begin(),neighbors.end());
    std::vector<kmer_id_t>::iterator it = std::unique_copy(
        neighbors.begin(), neighbors.end(), neighbors.begin());
    neighbors.resize(it - neighbors.begin());

#ifdef DEBUG
    std::cout << "Size : " << neighbors.size() << std::endl;
    for(int i = 0;i < neighbors.size();i++) {
        print_bits(neighbors.at(i),k*2);
    }
#endif

    free(a);
    return 0;
}

