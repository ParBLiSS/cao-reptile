#include "util.h"
#include "WorkDistribution.hpp"
#include <chrono>

struct WInt{
    int value;
    int sz;
    WInt(){
        value = 9;
        sz = 50;
        // std::cout << "C";
    }
    size_t size(){
        return sz;
    }
    void reset(){
        value = 0;
        sz = 0;
    }

    WInt(const WInt& other){
        //std::cout << other.value;
        value  = other.value;
        sz = other.sz;
    }
};

void swap(WInt& x, WInt& y){
    std::swap(x.value, y.value);
    std::swap(x.sz, y.sz);
}

struct WIntLoader{
    void operator()(const int& , int woffStart, int, WInt& tmp){
        tmp.value = woffStart;
    }
};

struct BatchIntWork{
    void operator()(const WInt& rbatch, int& ecr, int tid, int rank){
        int x = tid + rank * rbatch.value + ecr;
        x = 1;
        // sleep for quite a while
        std::this_thread::sleep_for(std::chrono::seconds(1));
        //std::cout << rank ;
    }
};

struct UnitIntWork{
    void operator()(const WInt& rbatch, int& ecr, unsigned ipos){
        int x = ipos * rbatch.value + ecr;
        x = x + 12;
        // sleep for a long while
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        //std::cout << "+" ;
    }
};

void check_file(const char* fname, int rank){
    std::ifstream fin(fname);
    if(!fin){
        if(rank == 0){
            std::cout << "ERROR:Can not open Input File!" << std::endl;
            std::cout << "Syntax:preptile /path/to/config-file" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

int main(int argc, char *argv[]){
    int provided, rs, rank, size;
    rs = MPI_Init_thread(&argc, &argv,
                         MPI_THREAD_FUNNELED, &provided );

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
        std::cout << (provided == MPI_THREAD_FUNNELED) << std::endl;

    //check_file(argv[1], rank);
    //auto tWork = getTotalWork(argv[1]);
    int nthreads = 2;
    int tWork = 5000;
    std::vector<int> pload;
    pload.resize(nthreads + 1);

    WorkDistribution<WInt, int, int, WIntLoader, BatchIntWork,
                     UnitIntWork> wds(tWork, pload, nthreads, 100);
    if(rank == 0) {
        wds.masterMain();
    } else {
        wds.slaveMain();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
        std::cout << std::endl;

    MPI_Finalize();
}
