#ifndef WORKDISTRIBUTION_H
#define WORKDISTRIBUTION_H

#include <mpi.h>
#include <iostream>
#include <thread>
#include <vector>
#include <queue>
#include <cassert>
#include <chrono>
#include <mutex>
#include <atomic>
#include <sstream>

#include "util.h"

template <typename T>
MPI_Datatype get_mpi_dt()
{
    throw std::runtime_error("Unsupported MPI datatype");
    // default to int
    return MPI_INT;
}

template <>
MPI_Datatype get_mpi_dt<double>(){
    return MPI_DOUBLE;
}
template <>
MPI_Datatype get_mpi_dt<int>(){
    return MPI_INT;
}

template <>
MPI_Datatype get_mpi_dt<unsigned>(){
    return MPI_UNSIGNED;
}

template <>
MPI_Datatype get_mpi_dt<unsigned long>(){
    return MPI_UNSIGNED_LONG;
}

template <typename T>
class SharedQueue{

public:
    SharedQueue() = default;
    SharedQueue(const SharedQueue&) = delete;
    SharedQueue& operator=(const SharedQueue&) = delete;
    bool pop(T& item){
        using std::swap;
        std::lock_guard<std::mutex> lock(m_qmut);
        if(m_queue.empty())
            return false;
        swap(m_queue.front(), item);
        m_queue.pop();
        return true;
    }

    void push(std::vector<T>& items){
        using std::swap;
        std::lock_guard<std::mutex> lock(m_qmut);
        for(auto ait = items.begin(); ait != items.end();++ait){
            T tmp;
            m_queue.push(tmp);
            swap(m_queue.back(), *ait);
        }
    }

    void push(T* items, int size){
        using std::swap;
        std::lock_guard<std::mutex> lock(m_qmut);
        for(int i = 0; i< size;i++){
            T tmp;
            m_queue.push(tmp);
            swap(m_queue.back(), items[i]);
        }
    }

    void push(T& item){
        using std::swap;
        std::lock_guard<std::mutex> lock(m_qmut);
        T tmp;
        m_queue.push(tmp);
        swap(m_queue.back(), item);
    }

    size_t size() const{
        return m_queue.size();
    }

    bool empty(){
        return m_queue.empty();
    }

private:
    std::queue<T> m_queue;
    std::mutex m_qmut;
};


enum WDState{
    ASSIGN_WORK = 1, // ready for work to be assigned
    PENDING_WORK,  // no more work left, but work queue not empty
    FINISHED_WORK // work queue empty, go and wait for thread for
};

enum WorkRequest{
    ASK_WORK_TAG = 0x1,
    SND_WORK_TAG = 0x2
};

template <typename WorkItemType, typename OffsetType, typename PayLoadType,
          typename BatchLoaderType, typename BatchWorkerType,
          typename UnitWorkerType>
class WorkDistribution{
    static SharedQueue<WorkItemType> wrkQueue;
    static std::atomic_bool wrkFinished;
    std::string fileName;
    OffsetType totalWork;
    OffsetType workChunk;
    unsigned numThreads;
    std::vector<PayLoadType>& payLoads;
    BatchLoaderType load_work;
    UnitWorkerType unit_work;
    int mpiSize;
    int mpiRank;

    //This function will be called from a thread,
    //   therfore any function that is called by this function should be thread safe
    static void worker_thread(int tid, int rank, PayLoadType& pload) {
        // my first two chunks
        //  rank * (2 * work_chunk) * numThreads;
        // std::cout << "Launched by thread\n";
        while(true) {
            WorkItemType work;
            if(wrkQueue.pop(work)){
                // TODO: do a batch of work
                // Right now, sleep the thread for a while
                std::stringstream out;
                //out << rank << " " << work << std::endl;
                std::this_thread::sleep_for(std::chrono::seconds(1));
                std::cout << rank ;
                BatchWorkerType batch_work;
                batch_work(tid, rank, work, pload);
            }
            if(wrkFinished.load(std::memory_order_relaxed))
                break;
        }
    }

    void startWorkers(std::vector<std::thread>& workers) {
        //Launch a group of worker threads
        for (unsigned i = 0; i < numThreads; ++i) {
            workers.push_back(std::thread(worker_thread, i, mpiRank,
                                          std::ref(payLoads[i + 1])));
        }
    }

    void joinWorkers(std::vector<std::thread>& workers) {
        // Wait for worker threads to finish it up
        for(auto wit = workers.begin(); wit != workers.end(); wit++){
            wit->join();
        }
    }

    void doWork(){
        // let coord thread to do some work
        static WorkItemType localWork;
        static unsigned workProgress = 0;
        // get the next work item from queue
        if(workProgress == 0 && !wrkQueue.pop(localWork)){
            return;
        }
        // TODO: do work of some predefined granularity
        unit_work(workProgress, localWork, payLoads[0]);
        workProgress += 1;
        if(workProgress >= localWork.size()){
            workProgress = 0;
            localWork.reset();
        }
    }

    // construct the work to be assigned as vector of offets
    void assignWork(std::vector<OffsetType>& wrkOffsets,
                     OffsetType& workAssigned){
        for(auto sit = wrkOffsets.begin(); sit != wrkOffsets.end(); sit++){
            if(workAssigned < totalWork){
                *sit = workAssigned;
                workAssigned += workChunk;
            } else {
                *sit = 0;
            }
        }
    }

    bool updateWorkQueue(std::vector<OffsetType>& wrkOffsets){
        // returns true when the whole of wrkOffsets is actual work!
        //  i.e. returns true when there is no 'zero work' assigned
        auto last = std::find(wrkOffsets.begin(), wrkOffsets.end(), 0);
        auto nwrk = std::distance(wrkOffsets.begin(), last);
        for(int i = 0; i < nwrk; i++){
            WorkItemType tmp;
            load_work(fileName, wrkOffsets[i],
                      wrkOffsets[i] + workChunk, tmp);
            wrkQueue.push(tmp);
            //std::cout << mpiRank << " " << nwrk << " "
            //<< wrkQueue.size() << std::endl;
        }

        return (nwrk == (long)wrkOffsets.size());
    }

    int probeQuery(MPI_Status& status){
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, ASK_WORK_TAG, MPI_COMM_WORLD, &flag, &status);
        return flag;
    }

    bool assignWorkRemote(OffsetType& wrkAssigned){
        // returns true when the whole of threadWork is actual work!
        //   or when nothing is assigned
        int count, recvMsg;
        MPI_Status status;
        MPI_Request request;
        std::vector<OffsetType> threadWork(numThreads);
        // if anybody asking work,
        //   either assign work to them or tell them there is no more work to do
        if(!probeQuery(status)){
            return true;
        }
        // who is asking for work ?
        MPI_Get_count(&status, MPI_INT, &count);
        assert(count == 1);
        MPI_Recv(&recvMsg, 1, MPI_INT, status.MPI_SOURCE,
                 ASK_WORK_TAG, MPI_COMM_WORLD, &status);
        // assign work proptional to numThreads
        assignWork(threadWork, wrkAssigned);
        // send the allocated
        MPI_Isend(&(threadWork[0]), numThreads, get_mpi_dt<OffsetType>(),
                  status.MPI_SOURCE, SND_WORK_TAG, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status); // should I need to wait ?

        return (threadWork.back() != 0);
    }

    void recvWorkRemote(std::vector<OffsetType>& recvMessage){
        // ask and recieve work from the master co-ord process
        MPI_Status status;
        MPI_Request request;
        int sendMsg = (int) numThreads;
        MPI_Isend(&sendMsg, 1, MPI_INT, 0, ASK_WORK_TAG, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        MPI_Recv(&(recvMessage[0]), numThreads, get_mpi_dt<OffsetType>(), 0,
                 SND_WORK_TAG, MPI_COMM_WORLD, &status);
    }

    void masterCoord(){
        // This function is not thread safe, only one thread in a process
        //  is allowed to run this function.
        static int nwrkZero = 1;
        static OffsetType wrkAssigned = mpiSize * workChunk;
        static WDState coordState = ASSIGN_WORK;

        switch(coordState){
        case ASSIGN_WORK:
            // Assign to work local threads, if they don't have enough
            if(wrkQueue.size() < 2 * numThreads){
                std::vector<OffsetType> threadWork(numThreads);
                assignWork(threadWork, wrkAssigned);
                if(!updateWorkQueue(threadWork)) // if assigned 'zero work'
                    coordState = PENDING_WORK; // update state
            }
            // Assign work to remote processes, if they are asking for it
            if(!assignWorkRemote(wrkAssigned)){
                // Update state, if I have assigned 'zero work'
                coordState = PENDING_WORK;
                nwrkZero += 1;
            }
            break;
        case PENDING_WORK:
            // I, root, haven't recieved message from every one:
            //   Assign 'zero work', to any remote process asking for work
            if(nwrkZero < mpiSize && !assignWorkRemote(wrkAssigned)) {
                nwrkZero += 1;
            }
            // I, root, have assigned 'zero work' to every process;
            //   Also, my queue is empty. I can finish work now.
            if(nwrkZero == mpiSize && wrkQueue.empty()){
                coordState = FINISHED_WORK;
            }
            break;
        case FINISHED_WORK:
            // may be do compare and exchange
            if(!wrkFinished.load(std::memory_order_relaxed)){
                wrkFinished.store(true, std::memory_order_relaxed);
            }
            break;
        }
    }

    int slaveCoord(){
        // This function is not thread safe, only one thread in a process
        //  is allowed to run this function.
        static WDState coordState = ASSIGN_WORK;
        switch(coordState){
        case ASSIGN_WORK:
            // If my local threads don't have enough work, ask from the root
            if(wrkQueue.size() < 2 * numThreads){
                std::vector<OffsetType> recvMessage(numThreads);
                recvWorkRemote(recvMessage);
                // If I have been assigned 'zero work', update state to pending
                if(!updateWorkQueue(recvMessage))
                    coordState = PENDING_WORK;
            }
            break;
        case PENDING_WORK:
            // No more work available, but work queue might not be empty
            if(wrkQueue.empty()){
                coordState = FINISHED_WORK;
            }
            break;
        case FINISHED_WORK:
            // No more work available and work queue is empty too
            if(!wrkFinished.load(std::memory_order_relaxed)){
                wrkFinished.store(true, std::memory_order_relaxed);
            }
            break;
        }
        return 0;
    }

public:
    void masterMain(){
        // This function is not thread safe, only one thread should be
        //  allowed to run this function.
        std::vector<std::thread> workers;
        wrkFinished.store(false, std::memory_order_relaxed);
        startWorkers(workers);
        do {
            doWork();
            masterCoord();
            if(wrkFinished.load(std::memory_order_relaxed))
                break;
        } while(true);
        joinWorkers(workers);
    }

    void slaveMain(){
        // This function is not thread safe, only one thread should be
        //  allowed to run this function.
        std::vector<std::thread> workers;
        wrkFinished.store(false, std::memory_order_relaxed);
        startWorkers(workers);
        do {
            doWork();
            slaveCoord();
            if(wrkFinished.load(std::memory_order_relaxed))
                break;
        } while(true);
        joinWorkers(workers);
    }

    WorkDistribution(OffsetType tWork, std::vector<PayLoadType>& refPayload,
                     unsigned nThreads = 2, OffsetType wChunk = 0):
        totalWork(tWork), numThreads(nThreads), payLoads(refPayload)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        if(wChunk > 0){
            workChunk = wChunk;
        } else {
            workChunk = totalWork/(mpiSize * numThreads);
        }
        assert(numThreads > 0);
        assert(workChunk > 0);
        assert(payLoads.size() == numThreads + 1);
    }
};

template <typename WorkItemType, typename OffsetType, typename PayLoadType,
          typename BatchLoaderType, typename BatchWorkerType,
          typename UnitWorkerType>
SharedQueue<WorkItemType> WorkDistribution<WorkItemType, OffsetType, PayLoadType,
                                           BatchLoaderType, BatchWorkerType,
                                           UnitWorkerType>::wrkQueue;

template <typename WorkItemType, typename OffsetType, typename PayLoadType,
          typename BatchLoaderType, typename BatchWorkerType,
          typename UnitWorkerType>
std::atomic_bool WorkDistribution<WorkItemType, OffsetType, PayLoadType,
                                           BatchLoaderType, BatchWorkerType,
                                           UnitWorkerType>::wrkFinished;


#endif // WORKDISTRIBUTION_H
