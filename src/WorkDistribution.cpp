//Create a group of C++11 threads from the main program

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

template <typename T>
class SharedQueue{
public:
    SharedQueue() = default;
    SharedQueue(const SharedQueue&) = delete;
    SharedQueue& operator=(const SharedQueue&) = delete;
    bool pop(T& item){
        std::lock_guard<std::mutex> lock(m_qmut);
        if(m_queue.empty())
            return false;
        item = m_queue.front();
        m_queue.pop();
        return true;
    }

    void push(const std::vector<T>& items){
        std::lock_guard<std::mutex> lock(m_qmut);
        for(auto ait = items.begin(); ait != items.end();++ait)
            m_queue.push(*ait);
    }

    void push(const T* items, int size){
        std::lock_guard<std::mutex> lock(m_qmut);
        for(int i = 0; i< size;i++)
            m_queue.push(items[i]);
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

enum WState{
    ASSIGN_WORK = 1, // ready for work to be assigned
    PENDING_WORK,  // no more work left, but work queue not empty
    FINISHED_WORK // work queue empty, go and wait for thread for
};

class WorkDistribution{
    static const int ASK_WORK_TAG = 200;
    static const int SND_WORK_TAG = 300;
    static const int total_work = 5000;
    static const int num_threads = 2;
    static const int work_chunk = 100;
    static SharedQueue<int> wrkQueue;
    static std::atomic_bool wrkFinished;
    //static WState work_state;
    int size;
    int rank;

    //This function will be called from a thread,
    //   therfore any function that is called by this function should be thread safe
    static void worker_thread(int rank) {
        // my first two chunks
        //  rank * (2 * work_chunk) * num_threads;
        // std::cout << "Launched by thread\n";
        while(true) {
            int work;
            if(wrkQueue.pop(work)){
                // do work
                // Right now, sleep the thread for a while
                std::stringstream out;
                out << rank << " " << work << std::endl;
                std::this_thread::sleep_for(std::chrono::seconds(1));
                std::cout << rank ;
            }
            if(wrkFinished.load(std::memory_order_relaxed))
                break;
        }
    }

    void start_workers(std::vector<std::thread>& workers) {
        //Launch a group of worker threads
        for (int i = 0; i < num_threads; ++i) {
            workers.push_back(std::thread(worker_thread, rank));
        }
    }

    void join_workers(std::vector<std::thread>& workers) {
        // Wait for worker threads to finish it up
        for(auto wit = workers.begin(); wit != workers.end(); wit++){
            wit->join();
        }
    }

    // construct the work to be assigned as vector of offets
    void assign_work(std::vector<int>& wrkOffsets, int& work_assigned){
        for(auto sit = wrkOffsets.begin(); sit != wrkOffsets.end(); sit++){
            if(work_assigned < total_work){
                *sit = work_assigned;
                work_assigned += work_chunk;
            } else {
                *sit = 0;
            }
        }
    }

    bool update_wqueue(std::vector<int>& wrkOffsets){
        // returns true when the whole of wrkOffsets is actual work!
        //  i.e. returns true when there is no 'zero work' assigned
        auto last = std::find(wrkOffsets.begin(), wrkOffsets.end(), 0);
        auto nwrk = std::distance(wrkOffsets.begin(), last);
        if(nwrk > 0){
            wrkQueue.push(&(wrkOffsets[0]), nwrk);
            //std::cout << rank << " " << nwrk << " "
            //<< wrkQueue.size() << std::endl;
        }

        return (nwrk == (long)wrkOffsets.size());
    }

    int probe_query(MPI_Status& status){
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, ASK_WORK_TAG, MPI_COMM_WORLD, &flag, &status);
        return flag;
    }

    bool assign_work_remote(int& wrkAssigned){
        // returns true when the whole of threadWork is actual work!
        //   or when nothing is assigned
        int count, recvMsg;
        MPI_Status status;
        MPI_Request request;
        std::vector<int> threadWork(num_threads);
        // if anybody asking work,
        //   either assign work to them or tell them there is no more work to do
        if(!probe_query(status)){
            return true;
        }
        // who is asking for work ?
        MPI_Get_count(&status, MPI_INT, &count);
        assert(count == 1);
        MPI_Recv(&recvMsg, 1, MPI_INT, status.MPI_SOURCE,
                 ASK_WORK_TAG, MPI_COMM_WORLD, &status);
        // assign work proptional to num_threads
        assign_work(threadWork, wrkAssigned);
        // send the allocated
        MPI_Isend(&(threadWork[0]), num_threads, MPI_INT, status.MPI_SOURCE,
                  SND_WORK_TAG, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status); // should I need to wait ?

        return (threadWork.back() != 0);
    }

    void recv_work_remote(std::vector<int>& recvMessage){
        MPI_Status status;
        MPI_Request request;
        int sendMsg = num_threads;
        MPI_Isend(&sendMsg, 1, MPI_INT, 0, ASK_WORK_TAG, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        MPI_Recv(&(recvMessage[0]), num_threads, MPI_INT, 0,
                 SND_WORK_TAG, MPI_COMM_WORLD, &status);
    }


    void send_work(){
        // This function is not thread safe, only one thread can call this fn.
        static int nwrkZero = 1;
        static int wrkAssigned = size * work_chunk;
        static WState distState = ASSIGN_WORK;
        bool updWork = true;

        switch(distState){
        case ASSIGN_WORK:
            // Assign to work local threads, if they don't have enough
            if(wrkQueue.size() < 2 * num_threads){
                std::vector<int> threadWork(num_threads);
                assign_work(threadWork, wrkAssigned);
                if(!update_wqueue(threadWork)) // if assigned 'zero work'
                    distState = PENDING_WORK; // update state
            }
            // Assign work to remote processes, if they are asking for it
            if(!assign_work_remote(wrkAssigned)){
                // Update state, if I have assigned 'zero work'
                distState = PENDING_WORK;
                nwrkZero += 1;
            }
            break;
        case PENDING_WORK:
            // I, root, haven't recieved message from every one:
            //   Assign 'zero work', to any remote process asking for work
            if(nwrkZero < size && !assign_work_remote(wrkAssigned)) {
                nwrkZero += 1;
            }
            // I, root, have assigned 'zero work' to every process;
            //   Also, my queue is empty. I can finish work now.
            if(nwrkZero == size && wrkQueue.empty()){
                distState = FINISHED_WORK;
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

    int recv_work(){
        // This function is not thread safe, only one thread can call this fn.
        static WState distState = ASSIGN_WORK;
        switch(distState){
        case ASSIGN_WORK:
            // If my local threads don't have enough work, ask from the root
            if(wrkQueue.size() < 2 * num_threads){
                std::vector<int> recvMessage(num_threads);
                recv_work_remote(recvMessage);
                // If I have been assigned 'zero work', update state to pending
                if(!update_wqueue(recvMessage))
                    distState = PENDING_WORK;
            }
            break;
        case PENDING_WORK:
            // No more work available, but work queue might not be empty
            if(wrkQueue.empty()){
                distState = FINISHED_WORK;
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
    void master_main(){
        // This function is not thread safe, only one thread can call this fn.
        std::vector<std::thread> workers;
        wrkFinished.store(false, std::memory_order_relaxed);
        start_workers(workers);
        do {
            send_work();
            if(wrkFinished.load(std::memory_order_relaxed))
                break;
        } while(true);
        join_workers(workers);
    }

    void slave_main(){
        // This function is not thread safe, only one thread can call this fn.
        std::vector<std::thread> workers;
        wrkFinished.store(false, std::memory_order_relaxed);
        start_workers(workers);
        do {
            recv_work();
            if(wrkFinished.load(std::memory_order_relaxed))
                break;
        } while(true);
        join_workers(workers);
    }

    WorkDistribution(){
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
};

SharedQueue<int> WorkDistribution::wrkQueue;
std::atomic_bool WorkDistribution::wrkFinished;

int main(int argc, char *argv[]){
    int provided, rs, rank, size;
    rs = MPI_Init_thread(&argc, &argv,
                         MPI_THREAD_FUNNELED, &provided );

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0)
        std::cout << provided << " " << MPI_THREAD_FUNNELED << std::endl;

    WorkDistribution wds;
    if(rank == 0) {
        wds.master_main();
    } else {
        wds.slave_main();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
        std::cout << std::endl;

    MPI_Finalize();
}
