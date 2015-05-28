#include "ECRunStats.hpp"
#include <sstream>
#include <fstream>

double elapsed(clock_t& end, clock_t& start){
  return (double (end - start))/ ((double) CLOCKS_PER_SEC);
}

double elapsed_local(timespec& finish, timespec& start){
  double tdiff;
  tdiff = (finish.tv_sec - start.tv_sec);
  tdiff += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
  return tdiff;
}


ECRunStats::ECRunStats(){
    tstartInit = MPI_Wtime();
    tstart = tstartInit;
}

void ECRunStats::reportTimings(Para& params, std::ostream& ofs){
    // Output for counting the number of failures and success
    //std::stringstream out;
    //out << params.oErrName << params.m_rank ;
    //ecdata.writeQueryStats(out.str());

    int p = params.m_size;
    if(params.m_rank == 0){
        std::stringstream oss;
        oss << "--" << std::endl
            << "nproc" << "\t" << "phase" << "\t"
            << "start" << "\t" << "stop" << "\t"
            << "duration" << std::endl;
        oss << p << "\t" << "read global" << "\t"
            << read_sync_start << "\t"
            << read_sync_stop << "\t"
            << read_sync_stop - read_sync_start
            << std::endl;
        oss << p << "\t" << "kmer global" << "\t"
            << kmer_sync_start << "\t"
            << kmer_sync_stop << "\t"
            << kmer_sync_stop - kmer_sync_start
            << std::endl;
        oss << p << "\t" << "ec global" << "\t"
            << ec_sync_start << "\t"
            << ec_sync_stop << "\t"
            << ec_sync_stop - ec_sync_start
            << std::endl;
        oss << p << "\t" << "final global" << "\t"
            << tstartInit << "\t" << tstop << "\t"
            << (tstop - tstartInit) << std::endl;
        oss << "--" << std::endl;
        oss << "rank" << "\t" << "phase" << "\t"
            << "start" << "\t" << "stop" << "\t"
            << "duration"
            << std::endl;
        ofs << oss.str();
        ofs.flush();
    }
    for(int i = 0; i < p; i++){
        MPI_Barrier(MPI_COMM_WORLD);
        if(i == params.m_rank){
            std::stringstream oss;
            oss << i << "\t" << "read local" << "\t"
                << elapsed_local(tstart_read_p, tstart_read_p) << "\t"
                << elapsed_local(tstop_read_p, tstart_read_p) << "\t"
                << elapsed_local(tstop_read_p, tstart_read_p)
                << std::endl;
            oss << i << "\t" << "kmer local" << "\t"
                << elapsed_local(tstart_kmer_p, tstart_read_p) << "\t"
                << elapsed_local(tstop_kmer_p, tstart_read_p) << "\t"
                << elapsed_local(tstop_kmer_p, tstart_kmer_p)
                << std::endl;
            oss << i << "\t" << "ec local" << "\t"
                << elapsed_local(tstart_ec_p, tstart_read_p) << "\t"
                << elapsed_local(tstop_ec_p, tstart_read_p) << "\t"
                << elapsed_local(tstop_ec_p, tstart_ec_p)
                << std::endl;
            ofs << oss.str();
            ofs.flush();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void ECRunStats::reportQueryCounts(ECData& ecdata, std::ostream& ofs){
    Para& params = ecdata.getParams();

    if(params.m_rank == 0) {
      std::stringstream oss;
      oss << "--" << std::endl;
      oss << "proc" << "\t" << "type" << "\t" << "query counts"  << "\t"
          << "query fails" << "\t" << "query success" << std::endl;
      ofs << oss.str();
      ofs.flush();
    }
    int p = params.m_size;
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i = 0; i < p; i++){
      MPI_Barrier(MPI_COMM_WORLD);
      if(i == params.m_rank){
        std::stringstream oss;
        oss << i << "\t" << "kmer" << "\t" << ecdata.getKmerQueries()
            << "\t" << ecdata.getKmerQueryFails() << "\t"
            << (ecdata.getKmerQueries()) - (ecdata.getKmerQueryFails())
            << std::endl;
        oss << i << "\t" << "tile" <<  "\t" << ecdata.getTileQueries()
            << "\t" << ecdata.getTileQueryFails() << "\t"
            << (ecdata.getTileQueries()) - (ecdata.getTileQueryFails())
            << std::endl;
        ofs << oss.str();
        ofs.flush();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(params.m_rank == 0) {
        std::stringstream oss;
        oss << "proc" << "\t" << "type" << "\t" ;
        for(unsigned j = 0; j < MAX_LEVELS; j++)
            oss << "L" << j << "\t" ;
        oss <<  std::endl;
        ofs << oss.str();
        ofs.flush();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for(int i = 0; i < p; i++){
      MPI_Barrier(MPI_COMM_WORLD);
      if(i == params.m_rank){
        std::stringstream oss;
        oss << i << "\t" << "kmer";
        for(unsigned j = 0; j < MAX_LEVELS; j++)
            oss << "\t" << ecdata.getKmerLevels(j);
        oss << std::endl;
        oss << i << "\t" << "tile";
        for(unsigned j = 0; j < MAX_LEVELS; j++)
            oss << "\t" << ecdata.getTileLevels(j);
        oss << std::endl;
        ofs << oss.str();
        ofs.flush();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void ECRunStats::updateFileReadTime(std::ostream&){
    read_sync_start = tstart;
    read_sync_stop = tstop;
}

void ECRunStats::updateDistSpectrumTime(ECData& ecdata, std::ostream& ofs){
    long kcount = ecdata.getKmerCount();
    long tilecount = ecdata.getTileCount();
    long tmp = kcount;
    MPI_Reduce( &tmp, &kcount, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
    tmp = tilecount;
    MPI_Reduce( &tmp, &tilecount, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );

    if(ecdata.getParams().m_rank != 0)
      return;

    std::stringstream oss, oss2;
    oss << "kmer count\t" << kcount << std::endl;
    oss << "tile count\t" << tilecount << std::endl;
    oss << "absent kmer\t" << ecdata.getParams().absentKmers << std::endl;
    std::cout << oss.str();
    std::cout.flush();

    kmer_sync_start = tstart;
    kmer_sync_stop = tstop;
    oss2 << "spectrum construction\t" << tstop-tstart << std::endl;
    ofs << oss2.str();
    ofs.flush();
}


void ECRunStats::updateSpectrumTime(ECData& ecdata, std::ostream& ofs){
    std::stringstream oss, oss2;
    oss << "kmer count\t" << ecdata.getKmerCount() << std::endl;
    oss << "tile count\t" << ecdata.getTileCount() << std::endl;
    oss << "absent kmer\t" << ecdata.getParams().absentKmers << std::endl;
    std::cout << oss.str();
    std::cout.flush();

    kmer_sync_start = tstart;
    kmer_sync_stop = tstop;
    oss2 << "spectrum construction\t" << tstop-tstart << std::endl;
    ofs << oss2.str();
    ofs.flush();
}

void ECRunStats::updateECTime(std::ostream& ofs){
    std::stringstream oss;
    ec_sync_start = tstart;
    ec_sync_stop = tstop;
    oss << "error correction\t" << tstop-tstart << std::endl;
    oss << "total\t" << tstop-tstartInit << std::endl;
    ofs << oss.str();
    ofs.flush();
}
