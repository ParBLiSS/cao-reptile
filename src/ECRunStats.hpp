#ifndef ECRUNSTATS_H
#define ECRUNSTATS_H

#include <fstream>
#include <time.h>

#include "util.h"
#include "ECData.hpp"

struct ECRunStats{
    double tstartInit, tstart, tstop;
    double read_sync_start, read_sync_stop,
        kmer_sync_start, kmer_sync_stop,
        ec_sync_start, ec_sync_stop;
    // clock_t tstart_read_p, tstop_read_p,
    //     tstart_kmer_p, tstop_kmer_p,
    //     tstart_ec_p, tstop_ec_p;
    timespec tstart_read_p, tstop_read_p,
        tstart_kmer_p, tstop_kmer_p,
        tstart_ec_p, tstop_ec_p;
    // update global timings
    void updateFileReadTime(std::ostream& ofs);
    void updateSpectrumTime(ECData& ecd, std::ostream& ofs);
    void updateECTime(std::ostream& ofs);
    void updateDistSpectrumTime(ECData& ecdata, std::ostream& ofs);
    void updateKmerDistSpectrum(ECData& ecdata, std::ostream& ofs);
    // report timigs
    void reportTimings(Para& params, std::ostream& ofs);
    void reportQueryCounts(ECData& ecd, std::ostream& ofs);
    void reportDynamicLoadTimings(ECData& ecd, std::vector<double>& stTimings,
                                  std::ostream& ofs);
    ECRunStats();
};

#endif /* ECRUNSTATS_H */
