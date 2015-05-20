#ifndef _KMER_COUNT_H
#define _KMER_COUNT_H

#include "ECData.hpp"

void count_kmers(ECData& ecdata);
void local_tile_spectrum(ECData& ecdata, bool qFlag = false);
void local_kmer_spectrum(ECData& ecdata, bool qFalg = true);

#endif
