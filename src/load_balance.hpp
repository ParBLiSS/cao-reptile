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

#ifndef _LOAD_BALANCE_H
#define _LOAD_BALANCE_H

#include <mpi.h>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <climits>
#include <sstream>

template <typename StructDataType>
void load_balance(StructDataType *&karray, int &kcount,int &ksize,
                  MPI_Datatype mpi_struct_type,
                  int size,int rank){
    int i, j;
    unsigned totalelements = 0;

    if(size == 1){
        return ;
    }
    int *pcount = new int[size](); // how much each processor has

    MPI_Allgather(&kcount, 1 , MPI_INT,
                  pcount, 1 , MPI_INT, MPI_COMM_WORLD);

    // Count the total number of elements
    for(i =0; i<size; i++){
        totalelements += pcount[i];
    }
    int tmp = totalelements/size;
    j = totalelements%size;

    int *recvcts =  new int[size]();
    int *sendcts = new int[size]();
    int *recvdisp = new int[size]();
    int *senddisp = new int[size]();
    int *underhood = new int[size]();
    int *overhood = new int[size]();
    int *elements_bloc = new int[size]();
    int *under = new int[size]();
    int *over = new int[size]();

    for( i =0; i<size; i++){
        if(i < j )
            elements_bloc[i] = tmp +1;
        else
            elements_bloc[i] = tmp;
        // Sriram : I DONT THINK THIS IS NEEDED!
        // if( rank == i)
        // newarray = (int *)malloc(elements_bloc[i]*sizeof(int));

        if(pcount[i] < elements_bloc[i]){
            under[i] = elements_bloc[i] - pcount[i];
            over[i] = 0;
        }
        else{
            under[i] = 0;
            over[i] = pcount[i] - elements_bloc[i];
        }
    }

#ifdef DEBUG
    std::stringstream out1;
    out1 << "PROC " << rank << " COUNT : " << pcount[rank]
         << " : UNDER COUNT : " << under[rank]
         << " : OVER COUNT : " << over[rank] << std::endl;
    std::cout << out1.str();
#endif

    underhood[0] = under[0];
    overhood[0] = over[0];

    for( i =1; i<size; i++){
        underhood[i] = underhood[i -1] + under[i];
        overhood[i] = overhood[i -1] + over[i];
    }

    tmp =0;
    if (under[rank] > 0){
        for(i=0; i<size; i++){
            sendcts[i]=0;
            senddisp[i]=0;
        }
        sendcts[rank] = kcount;
        j = underhood[rank] - under[rank] ;
        for(i=0;i<size;i++){
            if( i == rank){
                recvcts[i] = kcount;
                recvdisp[i] = tmp;
                tmp += kcount;
            }
            else{
                if(overhood[i] > j) {
                    if(tmp < elements_bloc[rank]){
                        if(rank < i){
                            if(tmp + overhood[i] -j <=elements_bloc[rank]){  // need to make sure that no one takes more than under[rank] elements. This depends on whether rank has appeared before i or will be appearing after i.
                                recvcts[i] = overhood[i] - j;
                                recvdisp[i]=tmp;
                                tmp += recvcts[i];
                                j = overhood[i];
                            }
                            else{
                                recvcts[i] = elements_bloc[rank] - tmp;
                                recvdisp[i] = tmp;
                                tmp = elements_bloc[rank];
                                j = overhood[i]; // does not matter whatever is set
                            }
                        } // close of if rank < i
                        else{
                            if(tmp + overhood[i] -j + kcount <= elements_bloc[rank]){
                                recvcts[i] = overhood[i] - j;
                                recvdisp[i]=tmp;
                                tmp += recvcts[i];
                                j = overhood[i];
                            }
                            else{
                                recvcts[i] = under[rank] - tmp;
                                recvdisp[i] = tmp;
                                tmp = elements_bloc[rank] - kcount;
                                j = overhood[i]; // need to check - ok . does not matter here as well
                            }
                        }  // close of else of if rank <i
                    }// close of if tmp < elements_bloc[rank]
                    else{
                        recvcts[i] = 0;
                        recvdisp[i] = 0;
                    }
                }// close of if overhood[i] > j
                else{
                    recvcts[i]  = 0;
                    recvdisp[i] = 0;
                }
            } //close of i != rank

        } //close of for loop

    }  // close of under[rank] > 0
    else{  //i.e under[rank] = 0
        for(i=0; i<size; i++){
            recvcts[i]=0;
            recvdisp[i]=0;
        }
        recvcts[rank] = elements_bloc[rank];
        j = overhood[rank] - over[rank] ;
        for(i=0;i<size;i++){
            if( i == rank){
                sendcts[i] = elements_bloc[rank];
                senddisp[i] = tmp;
                tmp += elements_bloc[rank];
            }
            else{
                if(underhood[i] > j) {
                    if(tmp < kcount){
                        if(rank < i){
                            if(tmp + underhood[i] -j <= kcount){   // need to make sure that no one takes more than under[rank] elements. This depends on whether rank has appeared before i or will be appearing after i.
                                sendcts[i] = underhood[i] - j;
                                senddisp[i]=tmp;
                                tmp += sendcts[i];
                                j = underhood[i];
                            }
                            else{
                                sendcts[i] = kcount - tmp;
                                senddisp[i] = tmp;
                                tmp = elements_bloc[rank];
                                j = underhood[size -1]; // does not matter whatever is set
                            }
                        } // close of if rank < i
                        else{
                            if(tmp + underhood[i] -j + elements_bloc[rank] <= kcount){
                                sendcts[i] = underhood[i] - j;
                                senddisp[i]=tmp;
                                tmp += sendcts[i];
                                j = underhood[i];
                            }
                            else{
                                sendcts[i] = over[rank] - tmp;
                                senddisp[i] = tmp;
                                tmp = over[rank];
                                j = underhood[size -1]; // need to check - ok . does not matter here as well
                            }
                        }  // close of else of if rank <i
                    }// close of if tmp < elements_bloc[rank]
                    else{
                        sendcts[i] = 0;
                        senddisp[i] = 0;
                    }
                }// close of if overhood[i] > j
                else{
                    sendcts[i]  = 0;
                    senddisp[i] = 0;
                }
            } //close of i != rank

        } //close of for loop

    }  // close of under[rank]  =  0

#ifdef DEBUG
    // just printing to test
    for(i = 0 ; i< size; i++){
        std::stringstream sout;
        if(rank == i){
            sout << "PROC : " << rank << std::endl
                 << "SNDCTS\tSNDDISP\tRCVCTS\tRCVDISP (BY COL)" << std::endl;
            for(j =0; j< size; j++)
                sout << sendcts[j] << "\t"
                     << senddisp[j] << "\t"
                     << recvcts[j] << "\t"
                     << recvdisp[j] << std::endl;
            sout << std::endl;
            std::cout << sout.str();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    StructDataType *newarray;
    newarray = (StructDataType *)malloc(elements_bloc[rank]*sizeof(StructDataType));

    MPI_Alltoallv(karray,sendcts,senddisp,mpi_struct_type,
                  newarray,recvcts,recvdisp,mpi_struct_type,MPI_COMM_WORLD);
#ifdef DEBUG
    for(i=0; i< size; i++){
        if(rank ==i){
            std::stringstream sout;
            sout << "PROC : " << rank << " TOTAL ELTS : " << elements_bloc[i] << std::endl;
            // for(j=0;j<elements_bloc[i];j++)
            //     sout << newarray[j].ID << "\t" ;
            sout << std::endl;
            std::cout << sout.str();
        }
        //MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
    free(karray);
    karray = newarray; ksize = kcount = elements_bloc[rank];

    delete[] recvcts;
    delete[] recvdisp;
    delete[] sendcts;
    delete[] senddisp;
    delete[] underhood;
    delete[] overhood;
    delete[] elements_bloc;
    delete[] under;
    delete[] over;
    delete[] pcount;
}

template <typename StructDataType>
void load_balance2(StructDataType *&karray, int &kcount,int &ksize,
                  MPI_Datatype mpi_struct_type,
                  int size,int rank){
    int i, j;
    unsigned totalelements = 0;
    if(size == 1){
        return ;
    }
    int *pcount = new int[size](); // how much each processor has

    MPI_Allgather(&kcount, 1 , MPI_INT,
                  pcount, 1 , MPI_INT, MPI_COMM_WORLD);

    bool lBalanceReqd = false;
    for (i=0;i<size;i++) {
        if(pcount[i] < size) {
            lBalanceReqd = true;
            break;
        }
    }
    if(!lBalanceReqd) {
        if (rank == 0) {
            std::cout << "load balancing\tNo" << std::endl;
        } else {
            std::cout << "load balancing\tYes " << std::endl;
        }
        std::cout.flush();
        return;
    }
    // Count the total number of elements
    for(i =0; i<size; i++){
        totalelements += pcount[i];
    }

    int *recvcts =  new int[size]();
    int *sendcts = new int[size]();
    int *recvdisp = new int[size]();
    int *senddisp = new int[size]();
    int *elements_bloc = new int[size]();
    int *under = new int[size]();
    int *over = new int[size]();

    //---------------------------------------------------------
    int *sendproc = new int[size]();  // ranks of proc that send
    int *recvproc = new int[size]();  // ranks of proc that recv
    int sndproccount = 0, // number of processes that only send i.e. are over
        recvproccount = 0; // number of processes that only rcv i.e. are under
    int tmp = totalelements/size;
    j = totalelements%size;
    // Compute the over-under for each process
    for( i =0; i<size; i++){
        if(i < j )
            elements_bloc[i] = tmp +1;
        else
            elements_bloc[i] = tmp;

        if(pcount[i] < elements_bloc[i]){
            under[i] = elements_bloc[i] - pcount[i];
            over[i] = 0;
            recvproc[recvproccount] = i;
            recvproccount++;
        }
        else{
            under[i] = 0;
            over[i] = pcount[i] - elements_bloc[i];
            sendproc[sndproccount] = i;
            sndproccount++;
        }
    }

    int  my_sent_offset = 0, // offset starting from which i should send
        my_rcvd_offset = 0; // offset starting from which i should rcv
    // Set up the number of elements I should recv from myself
    if(over[rank] > under[rank]){
        recvcts[rank] = sendcts[rank] = elements_bloc[rank];
        recvdisp[rank] = senddisp[rank] = 0;
        my_sent_offset = elements_bloc[rank];
    } else {
        recvcts[rank] = sendcts[rank] = pcount[rank];
        recvdisp[rank] = senddisp[rank] = 0;
        my_rcvd_offset = pcount[rank];
    }

    // With a merge-like algo. set-up the send/recv counts
    int sndptr = 0, rcvptr = 0, sent = 0,recvd = 0;
    while(sndptr < sndproccount && rcvptr < recvproccount) {
        int cursndrank = sendproc[sndptr],
            currcvrank = recvproc[rcvptr];
        int currsend_avail = over[cursndrank] - sent,
            currecv_avail = under[currcvrank] - recvd;
        int xfernow  = 0;
        // xfer the min of currecv_avail, currsend_avail
        if(currecv_avail <= currsend_avail) {
            xfernow = currecv_avail;
        }
        else{
            xfernow  = currsend_avail;
        }
        // If I am going to send, then set up the send counts
        //  against the recv rank process
        if(rank == cursndrank) {
            sendcts[currcvrank] = xfernow;
            senddisp[currcvrank] = my_sent_offset;
            my_sent_offset += xfernow;
        }
        // If I am going to recv, then set up the rcv counts
        //  against the send rank process
        if(rank == currcvrank) {
            recvcts[cursndrank] = xfernow;
            recvdisp[cursndrank] = my_rcvd_offset;
            my_rcvd_offset += xfernow;
        }
        // Xfer set-up done, now update snd,rcv pts and
        //  reset sent,recvd if necessary
        sent += xfernow;
        recvd += xfernow;
        if(sent == over[cursndrank]) {
            sndptr++;
            sent = 0;
        }
        if(recvd == under[currcvrank]) {
            rcvptr++;
            recvd = 0;
        }
    }
#ifdef DEBUG
    // just printing to test
    for(i = 0 ; i< size; i++){
        std::stringstream sout;
        if(rank == i){
            sout << "PROCX : " << rank << std::endl
                 << "SNDCTS\tSNDDISP\tRCVCTS\tRCVDISP (BY COLX)" << std::endl;
            for(j =0; j< size; j++)
                sout << sendcts[j] << "\t"
                     << senddisp[j] << "\t"
                     << recvcts[j] << "\t"
                     << recvdisp[j] << std::endl;
            sout << std::endl;
            std::cout << sout.str();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
    //---------------------------------------------------------

    StructDataType *newarray;
    newarray = (StructDataType *)malloc(elements_bloc[rank]*sizeof(StructDataType));

    MPI_Alltoallv(karray,sendcts,senddisp,mpi_struct_type,
                  newarray,recvcts,recvdisp,mpi_struct_type,
                  MPI_COMM_WORLD);
    free(karray);
    karray = newarray; ksize = kcount = elements_bloc[rank];

    delete[] recvcts;
    delete[] recvdisp;
    delete[] sendcts;
    delete[] senddisp;
    delete[] elements_bloc;
    delete[] under;
    delete[] over;
    delete[] pcount;
    delete[] recvproc;
    delete[] sendproc;
}

#endif
