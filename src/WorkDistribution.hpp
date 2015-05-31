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

#include "mpi_util.hpp"
//
// A share queue using a lock
//  - T should have constructor with empty arguments
//  - swap(T&, T&) should be implemented to avoid unnecessary copying
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
    ASSIGN_WORK = 0, // ready for work to be assigned
    PENDING_WORK,  // no more work left, but work queue not empty
    FINISHED_WORK, // work queue empty, go and wait for thread for
    NR_WORK_STATES // place holder for counting states
};

enum WorkRequest{
    ASK_WORK_TAG = 0x1,
    SND_WORK_TAG = 0x2
};

// WorkDistribution
//
// - Distributes a given quantity of work into workchunks across the
//   all distributed processes and their shared memory threads.
// - Architected as a single master procecss and a bunch of slave
//   processes
// - Each process has a co-ordination thread and  'numWorkers' no. of
///  worker threads.
// - Initial work chunk is assigned without explicit allocation by the
//   master process.
// - Co-ordination threads in slave processes have the responsibility
//   to allocate work for worker threads by asking work from the
//   master process.
// - Co-ordination threads in master process have the responsibility
//   to allocate work for its worker threads and also allocate work for
//   all slave processes.
// - Co-ordination threads in all the processes also load the work item
//   (incase it needs to be loaded from a file or network) and do some
//   part of the work in between other responsibilites.
//
// WorkItemType holds the information regarding work that needs to be
// done and is expected to have the following interface
//   0. Constructor with no arguments
//   1. Has a size() function that returns
//   2. Has a reset() function that clears the state of the work item
//   3. implement a swap function to avoid copying in the queue
// size() should correspond to the number of units of work.
//
// SizeType is used to measure the quantity of work, size of the work
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
// SizeType startOffset, SizeType endOffset, WorkItemType& item
//
// BatchWokerType does the work corresponding to a given batch. Its
// arguments are WorkItemType& item, PayLoadType& pl, int threadid, int rank
//
// UnitWokerType does the work corresponding to a single unit with in a
// batch. Its arguments are WorkItemType& item, PayLoadType& pl,
// unsigned unit_id
//
template <typename WorkItemType, typename SizeType, typename PayLoadType,
          typename BatchLoaderType, typename BatchWorkerType,
          typename UnitWorkerType, typename ChunkSizeType,
          typename ParamType>
class WorkDistribution{
    SharedQueue<WorkItemType> wrkQueue;
    std::atomic_bool wrkFinished;

    // work distribution
    SizeType totalWork;
    std::vector<PayLoadType>& payLoads;
    unsigned numWorkers;
    SizeType initChunk; // initial allocation chunk size

    // woker functions
    BatchLoaderType load_work;
    UnitWorkerType unit_work;
    ChunkSizeType chunk_size;

    // work item and progress indicator for co-ordination thread
    WorkItemType crdWork;
    unsigned workProgress;

    // mpi information
    int mpiSize;
    int mpiRank;

    // co-ordination variables
    int numWorkZero; // master co-ord for counting no. 'zero work'
    SizeType workAssigned; // master co-ord tracks how much work is assinged
    WDState coordState; // current state

    std::vector<timespec> stateTimings;
    std::stringstream msgOut;

    const ParamType& inParams;

    //This function will be called from a thread,
    //   therfore any function that is called by this function should be thread safe
    // It is assumed that the PayLoadType object is exclusive to this thread
    static void worker_thread(int tid, int rank, PayLoadType& pload,
                              SharedQueue<WorkItemType>& wrkQueueRef,
                              std::atomic_bool& wrkFinishedRef) {
        BatchWorkerType batch_work;
        WorkItemType work;

        while(true) {
            if(wrkQueueRef.pop(work)){
                // do a batch of work
                batch_work(work, pload, tid, rank);
                work.reset();
            }
            if(wrkFinishedRef.load(std::memory_order_relaxed))
                break;
        }
    }

    // Launch a group of worker threads
    void startWorkers(std::vector<std::thread>& workers) {
        for (unsigned tid = 0; tid < numWorkers; ++tid) {
            // initialize the first chunk
            BatchLoaderType load_work;
            WorkItemType work;
            // My first work chunk
            //   (work_chunk) ((rank * (numWorkers + 1)) + tid);
            SizeType threadOffset = mpiRank * (numWorkers + 1) * initChunk;
            threadOffset += (tid + 1) * initChunk;

            bool initLoad = load_work(payLoads[tid + 1], threadOffset,
                                      threadOffset + initChunk, work);
            if(!initLoad)
                msgOut << "E";
            wrkQueue.push(work);
            // start thread
            workers.push_back(std::thread(worker_thread, tid, mpiRank,
                                          std::ref(payLoads[tid + 1]),
                                          std::ref(wrkQueue),
                                          std::ref(wrkFinished)));
        }
    }

    // Wait for worker threads to finish
    void joinWorkers(std::vector<std::thread>& workers) {
        for(auto wit = workers.begin(); wit != workers.end(); wit++){
            wit->join();
        }
    }

    // Work done by coord thread
    bool doWork(){
        // get the next work item from queue
        if(workProgress == 0 &&
           crdWork.size() == 0 && !wrkQueue.pop(crdWork)){
            return false;
        }
        // TODO: do work of some predefined granularity
        if(workProgress < crdWork.size()){
            unit_work(crdWork, payLoads[0], workProgress);
            workProgress += 1;
        } else {
            workProgress = 0;
            crdWork.reset();
        }
        return true;
    }

    // Compute the work to be assigned as vector of offets
    void assignWork(std::vector<SizeType>& wrkOffsets){
        for(auto sit = wrkOffsets.begin(); sit != wrkOffsets.end(); sit++){
            if(workAssigned < totalWork){
                *sit = workAssigned;
                workAssigned += chunk_size(totalWork, workAssigned, inParams);
            } else {
                *sit = 0;
            }
        }
    }

    // Update work queue with new set of items
    //   - returns true when the whole of wrkOffsets is actual work!
    //     i.e. returns true when there is no 'zero work' assigned
    bool updateWorkQueue(std::vector<SizeType>& wrkOffsets){
        auto last = std::find(wrkOffsets.begin(), wrkOffsets.end(), 0);

        for(auto ait = wrkOffsets.begin(); ait != last; ait++){
            WorkItemType tmp;
            if(load_work(payLoads[0], *ait,
                         (*ait) + chunk_size(totalWork, *ait, inParams), tmp))
                wrkQueue.push(tmp);
        }
        //std::cout << wrkOffsets.back();

        return (wrkOffsets.back() != 0);
    }

    // Check if any one is asking for work assignment
    int probeQuery(MPI_Status& status){
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, ASK_WORK_TAG, MPI_COMM_WORLD, &flag, &status);
        return flag;
    }

    // Assigns work to remote processes, if they are asking for it
    //  - returns true when the whole of threadWork is actual work!
    //    or when nothing is assigned
    bool assignWorkRemote(){
        int count, recvMsg, nSent = 0;
        MPI_Status status;
        MPI_Request request;
        std::vector<SizeType> threadWork(numWorkers);
        bool fullAssigned = true;
        // if anybody asking work,
        //   Either assign work to them or tell them there is no more work to do
        while(probeQuery(status)){
            // who is asking for work ?
            MPI_Get_count(&status, MPI_INT, &count);
            assert(count == 1);
            MPI_Recv(&recvMsg, 1, MPI_INT, status.MPI_SOURCE,
                     ASK_WORK_TAG, MPI_COMM_WORLD, &status);
            // always assign work of size 'numWorkers'
            assignWork(threadWork);
            // send the allocated work
            MPI_Isend(&(threadWork[0]), numWorkers, get_mpi_dt<SizeType>(),
                      status.MPI_SOURCE, SND_WORK_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status); // should I need to wait ?
            nSent += 1;
            if(threadWork.back() == 0){
                fullAssigned = false;
                numWorkZero += 1;
            }
        }

        return fullAssigned;
    }

    // Ask and recieve work from the master co-ord process
    void recvWorkRemote(std::vector<SizeType>& recvMessage){
        MPI_Status status;
        MPI_Request request;
        int sendMsg = (int) numWorkers;
        recvMessage.resize(numWorkers);
        MPI_Isend(&sendMsg, 1, MPI_INT, 0, ASK_WORK_TAG, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        MPI_Recv(&(recvMessage[0]), numWorkers, get_mpi_dt<SizeType>(), 0,
                 SND_WORK_TAG, MPI_COMM_WORLD, &status);
    }

    void updateCoordState(WDState newState){
        coordState = newState;
        if(coordState >= stateTimings.size())
            stateTimings.resize(coordState + 1);
        stateTimings[coordState] = local_time();
    }

    // Co-ordination thread in master process.
    //  - Fullfills the responsibilites based on the coordState
    //  - This function is not thread safe, only one thread in a process
    //     is allowed to run this function.
    void masterCoord(){

        switch(coordState){
        case ASSIGN_WORK:
            // Assign to work local threads, if they don't have enough
            if(wrkQueue.size() < 2 * numWorkers){
                std::vector<SizeType> threadWork(numWorkers);
                assignWork(threadWork);
                if(!updateWorkQueue(threadWork)) // if assigned 'zero work'
                    updateCoordState(PENDING_WORK); // update state
            }
            // Assign work to remote processes, if they are asking for it
            if(!assignWorkRemote()){
                // Update state, if I have assigned 'zero work'
                updateCoordState(PENDING_WORK);
            }
            break;
        case PENDING_WORK:
            // I, root, haven't recieved message from every one:
            //   Assign 'zero work', to any remote process asking for work
            if(numWorkZero < mpiSize){
                assignWorkRemote();
            }
            // I, root, have assigned 'zero work' to every process;
            //   Also, my queue is empty. I can finish work now.
            if(numWorkZero == mpiSize && wrkQueue.empty()){
                updateCoordState(FINISHED_WORK);
            }
            break;
        case FINISHED_WORK: // nothing to done, signal workers
            if(!wrkFinished.load(std::memory_order_relaxed)){
                wrkFinished.store(true, std::memory_order_relaxed);
            }
            break;
        default:
            break;
        }
    }

    // Co-ordination thread in slave process.
    //  - Fullfills the responsibilites based on the coordState
    //  - This function is not thread safe, only one thread in a process
    //     is allowed to run this function.
    int slaveCoord(){

        switch(coordState){
        case ASSIGN_WORK:
            // If my local threads don't have enough work, ask from the root
            if(wrkQueue.size() < 2 * numWorkers){
                std::vector<SizeType> recvMessage(numWorkers);
                recvWorkRemote(recvMessage);
                // If I have been assigned 'zero work', update state to pending
                if(!updateWorkQueue(recvMessage))
                    updateCoordState(PENDING_WORK);
            }
            break;
        case PENDING_WORK:
            // No more work available, but work queue might not be empty
            if(wrkQueue.empty()){
                updateCoordState(FINISHED_WORK);
            }
            break;
        case FINISHED_WORK:
            // No more work available and work queue is empty too
            if(!wrkFinished.load(std::memory_order_relaxed)){
                wrkFinished.store(true, std::memory_order_relaxed);
            }
            break;
        default:
            break;
        }
        return 0;
    }

    // Load a batch as the work item for coord thread
    void loadCoordWork(){
        SizeType startOffset = mpiRank * (numWorkers + 1) * initChunk;
        // std::cout << mpiRank << " " << startOffset << std::endl;
        bool cwork = load_work(payLoads[0], startOffset,
                               startOffset + initChunk, crdWork);
        if(!cwork)
            msgOut << "E";
    }

    // Inital work allocation
    void initQueue(){
        std::vector<SizeType> iOffsets(numWorkers + 1);
        SizeType currOffset;
        for(int j = 1; j <= INIT_QUEUE_FACTOR; j++){
            currOffset = mpiRank + (j * mpiSize);
            currOffset *= (numWorkers + 1) * initChunk;
            for(auto oit = iOffsets.begin(); oit != iOffsets.end(); ++oit){
                *oit = currOffset;
                currOffset += initChunk;
                // std::cout << mpiRank << " " << currOffset << std::endl;
            }
            updateWorkQueue(iOffsets);
        }
    }

public:

    const int INIT_QUEUE_FACTOR = 2;
    const std::vector<timespec>& getStateTimings() const{
        return stateTimings;
    }

    // Main loop of the co-rodnation thread
    //  - assumes the root process
    //  - this function is not thread safe, only the main thread should be
    //     allowed to run this function.
    void main(int root = 0){
        stateTimings.resize(NR_WORK_STATES + 1);
        wrkFinished.store(false, std::memory_order_relaxed);
        updateCoordState(ASSIGN_WORK);
        std::vector<std::thread> workers;
        startWorkers(workers); // start workers

        loadCoordWork(); // Load the local work first
        initQueue(); // Initializes the queue with 2 * numWorkers

        workAssigned = (1 + INIT_QUEUE_FACTOR) * mpiSize *
            (numWorkers + 1) * initChunk;
        // std::cout << mpiRank << " " << workAssigned << std::endl;
        // Start co-ordination loop
        do {
            doWork();
            if(mpiRank == root)
                masterCoord();
            else
                slaveCoord();
            if(wrkFinished.load(std::memory_order_relaxed))
                break;
        } while(true);

        while(doWork());  // finish if any local work left to do

        joinWorkers(workers);
        updateCoordState(NR_WORK_STATES);
        if(msgOut.str().size() > 0){
            msgOut << std::endl;
            std::cout << msgOut.str();
        }
    }

    WorkDistribution(SizeType tWork, std::vector<PayLoadType>& refPayload,
                     unsigned nWorkers, const ParamType& params):
        totalWork(tWork), payLoads(refPayload),
        numWorkers(nWorkers), inParams(params)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

        assert(numWorkers > 0);
        assert(payLoads.size() == numWorkers + 1);

        workProgress = 0;
        numWorkZero = 1;
        initChunk = chunk_size(totalWork, 0, inParams);
        assert(initChunk > 0);
    }
};

#endif // WORKDISTRIBUTION_H
