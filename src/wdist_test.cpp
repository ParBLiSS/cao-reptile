#include "util.h"
#include "WorkDistribution.hpp"

struct MyInt{
    int value;
    size_t size(){
        return 2;
    }
    void reset(){
        value = 0;
    }
};

struct MyIntLoader{
    void operator()(const std::string& fname, unsigned long woffStart,
                    unsigned long woffEnd, MyInt& tmp){
        tmp.value = woffStart;
    }
};

struct BatchIntWork{
    void operator()(int tid, int rank, const MyInt& rbatch, int& ecr){
        int x = tid + rank * rbatch.value + ecr;
        x = 1;
        // TODO: sleep ?
    }
};

struct UnitIntWork{
    void operator()(unsigned ipos, const MyInt& rbatch, int& ecr){
        int x = ipos * rbatch.value + ecr;
        x = x + 12;
        // TODO: sleep
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

unsigned long getTotalWork(const std::string& fileName){
    std::ifstream fin(fileName.c_str());
    fin.seekg(0,std::ios::end);
    return fin.tellg();
}

int main(int argc, char *argv[]){
    int provided, rs, rank, size;
    rs = MPI_Init_thread(&argc, &argv,
                         MPI_THREAD_FUNNELED, &provided );

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
        std::cout << provided << " " << MPI_THREAD_FUNNELED << std::endl;


    check_file(argv[1], rank);

    std::vector<int> pload;
    pload.resize(3);

    auto tWork = getTotalWork(argv[1]);

    WorkDistribution<MyInt, unsigned long, int,
                     MyIntLoader, BatchIntWork, UnitIntWork> wds(tWork, pload);
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
