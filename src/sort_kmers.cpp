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

#include "sort_kmers.hpp"
#include "ECData.hpp"
#include "util.h"


void kmer_sort(ECData *ecdata) {

    Para *params = ecdata->m_params;
    // kmer_t *allKmers;
    // int allKmerCount;
    // sort the k-mer pairs
    sort_kmers<kmer_t,kmer_id_t,KmerComp,KmerIDComp>(
        ecdata->m_karray,ecdata->m_kcount,ecdata->m_ksize,
        ecdata->m_mpi_kmer_t,mpi_kmer_id_t,
        KmerComp(),KmerIDComp(),
        true, params->tCard,
        ecdata->m_params);

    // tile_t *allTiles;
    // int allTileCount;
    // sort the tiles
    sort_kmers< tile_t,tile_id_t, TileComp, TileIDComp >(
        ecdata->m_tilearray,ecdata->m_tilecount, ecdata->m_tilesize,
        ecdata->m_mpi_tile_t,mpi_tile_id_t,
        TileComp(),TileIDComp(),
        true,params->tCard,
        ecdata->m_params);

    //
      ecdata->buildCacheOptimizedLayout();
}
