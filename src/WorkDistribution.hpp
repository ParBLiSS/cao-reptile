#ifndef WORKDISTRIBUTION_H
#define WORKDISTRIBUTION_H

#include <iostream>
#include <thread>
#include <vector>
#include <queue>
#include <cassert>
#include <mutex>
#include <atomic>
#include <sstream>
#include <mpi.h>

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
//
//
// WorkItemType holds the information regarding work that needs to be
// done and is expected to have the following interface
//   0. Constructor with no arguments
//   1. Has a size() function that returns
//   2. Has a reset() function that clears the state of the work item
//   3. implement a swap function to avoid copying in the queue
// size() should correspond to the number of units of work.
//
// OffsetType is used to measure the quantity of work, size of the work
// chunk, and offsets work. This type is expected to an integral type.
//
// A reference to PayLoadType object passsed to the worker
// functions. This object is expected to have all the information
// necessary to be made use of during the work. Each thread is supposed
// to get its own copy of a PayLoadType object so that the thread can
// write any results to this object with out any conflict.
//
// BatchLoaderType, BatchWorkerType and UnitWorkerType are function
// objects.
//
// BatchLoaderType loads a WorkItemType object from starting and ending
// offset. It also takes the payload object as an argument to retrieve
// the meta data necessary. The arguments are const PayLoadType& pl,
// OffsetType startOffset, OffsetType endOffset, WorkItemType& item
//
// BatchWokerType does the work corresponding to a given batch. Its
// arguments are WorkItemType& item, PayLoadType& pl, int threadid, int rank
//
// UnitWokerType does the work corresponding to a single unit with in a
// batch. Its arguments are WorkItemType& item, PayLoadType& pl,
// unsigned unit_id
//
template <typename WorkItemType, typename OffsetType, typename PayLoadType,
          typename BatchLoaderType, typename BatchWorkerType,
          typename UnitWorkerType>
class WorkDistribution{
    static SharedQueue<WorkItemType> wrkQueue;
    static std::atomic_bool wrkFinished;

    // work distribution
    OffsetType totalWork;
    std::vector<PayLoadType>& payLoads;
    unsigned numThreads;
    OffsetType workChunk;

    // woker functions
    BatchLoaderType load_work;
    UnitWorkerType unit_work;

    // work item and progress indicator for co-ordination thread
    WorkItemType crdWork;
    unsigned workProgress;

    // mpi information
    int mpiSize;
    int mpiRank;

    //This function will be called from a thread,
    //   therfore any function that is called by this function should be thread safe
    // It is assumed that the PayLoadType object is exclusive to this thread
    static void worker_thread(int tid, int nthreads, int rank,
                              OffsetType wChunk, PayLoadType& pload) {
        BatchLoaderType load_work;
        BatchWorkerType batch_work;
        WorkItemType work;
        // my first chunk
        //   (work_chunk) ((rank * (numThreads + 1)) + tid);
        OffsetType threadOffset = rank * (nthreads + 1) * wChunk;
        threadOffset += tid * wChunk;
        load_work(pload, threadOffset, threadOffset + wChunk, work);
        batch_work(work, pload, tid, rank);

        while(true) {
            work.reset();
            if(wrkQueue.pop(work)){
                // do a batch of work
                batch_work(work, pload, tid, rank);
            }
            if(wrkFinished.load(std::memory_order_relaxed))
                break;
        }
    }

    // Launch a group of worker threads
    void startWorkers(std::vector<std::thread>& workers) {
        for (unsigned i = 0; i < numThreads; ++i) {
            workers.push_back(std::thread(worker_thread, i, numThreads,
                                          mpiRank, workChunk,
                                          std::ref(payLoads[i + 1])));
        }
    }

    // Wait for worker threads to finish it up
    void joinWorkers(std::vector<std::thread>& workers) {
        for(auto wit = workers.begin(); wit != workers.end(); wit++){
            wit->join();
        }
    }

    // Work done by coord thread
    bool doWork(){
        // TODO: think about loading from work available
        // get the next work item from queue
        if(workProgress == 0 &&
           crdWork.size() == 0 && !wrkQueue.pop(crdWork)){
            return false;
        }
        // TODO: do work of some predefined granularity
        unit_work(crdWork, payLoads[0], workProgress);
        workProgress += 1;
        if(workProgress >= crdWork.size()){
            workProgress = 0;
            crdWork.reset();
        }
        return true;
    }

    // Compute the work to be assigned as vector of offets
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

    // Update work queue with new set of items
    //   - returns true when the whole of wrkOffsets is actual work!
    //     i.e. returns true when there is no 'zero work' assigned
    bool updateWorkQueue(std::vector<OffsetType>& wrkOffsets){
        auto last = std::find(wrkOffsets.begin(), wrkOffsets.end(), 0);
        //auto nwrk = std::distance(wrkOffsets.begin(), last);
        //for(int i = 0; i < nwrk; i++){
        for(auto ait = wrkOffsets.begin(); ait != last; ait++){
            WorkItemType tmp;
            load_work(payLoads[0], *ait, (*ait) + workChunk, tmp);
            wrkQueue.push(tmp);
        }

        //return (nwrk == (long)wrkOffsets.size());
        return (wrkOffsets.back() != 0);
    }

    int probeQuery(MPI_Status& status){
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, ASK_WORK_TAG, MPI_COMM_WORLD, &flag, &status);
        return flag;
    }

    // Assigns work to remote processes, if they are asking for it
    //  - returns true when the whole of threadWork is actual work!
    //    or when nothing is assigned
    bool assignWorkRemote(OffsetType& wrkAssigned){
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

    // Ask and recieve work from the master co-ord process
    void recvWorkRemote(std::vector<OffsetType>& recvMessage){
        MPI_Status status;
        MPI_Request request;
        int sendMsg = (int) numThreads;
        MPI_Isend(&sendMsg, 1, MPI_INT, 0, ASK_WORK_TAG, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        MPI_Recv(&(recvMessage[0]), numThreads, get_mpi_dt<OffsetType>(), 0,
                 SND_WORK_TAG, MPI_COMM_WORLD, &status);
    }

    // Co-ordination thread in master process.
    //  - Fullfills the responsibilites based on the coordState
    //  - This function is not thread safe, only one thread in a process
    //     is allowed to run this function.
    void masterCoord(){
        int nwrkZero = 1;
        OffsetType wrkAssigned = mpiSize * (numThreads + 1) * workChunk;
        WDState coordState = ASSIGN_WORK;

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
            // TODO: do compare and exchange
            if(!wrkFinished.load(std::memory_order_relaxed)){
                wrkFinished.store(true, std::memory_order_relaxed);
            }
            break;
        }
    }

    // Co-ordination thread in slave process.
    //  - Fullfills the responsibilites based on the coordState
    //  - This function is not thread safe, only one thread in a process
    //     is allowed to run this function.
    int slaveCoord(){
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

    // Load a batch as the work item for coord thread
    void loadWork(){
        OffsetType startOffset = mpiRank * (numThreads + 1) * workChunk;
        load_work(payLoads[0], startOffset,
                  startOffset + workChunk, crdWork);
    }

public:
    // Main loop of master co-rodnation thread
    //  - This function is not thread safe, only one thread should be
    //    allowed to run this function.
    void masterMain(){
        wrkFinished.store(false, std::memory_order_relaxed);
        std::vector<std::thread> workers;
        startWorkers(workers); // start workers

        loadWork(); // Load the local work first

        // Start co-ordination loop
        do {
            doWork();
            masterCoord();
            if(wrkFinished.load(std::memory_order_relaxed))
                break;
        } while(true);

        while(doWork());  // finish if any local work left to do

        joinWorkers(workers);
    }

    // Main loop of master co-rodnation thread
    //  - this function is not thread safe, only one thread should be
    //     allowed to run this function.
    void slaveMain(){
        wrkFinished.store(false, std::memory_order_relaxed);
        std::vector<std::thread> workers;
        startWorkers(workers); // start workers

        loadWork(); // Load the local work first

        do {
            doWork();
            slaveCoord();
            if(wrkFinished.load(std::memory_order_relaxed))
                break;
        } while(true);

        while(doWork());  // finish if any local work left to do

        joinWorkers(workers);
    }

    WorkDistribution(OffsetType tWork, std::vector<PayLoadType>& refPayload,
                     unsigned nThreads = 2, OffsetType wChunk = 0):
        totalWork(tWork), payLoads(refPayload),
        numThreads(nThreads), workChunk(wChunk)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        if(!(wChunk > 0)) {
            workChunk = totalWork/(mpiSize * numThreads);
        }
        workProgress = 0;
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
