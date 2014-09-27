#include <vector>
#include <iostream>
#include <iterator>
#include <cassert>
#include <algorithm>
#include <stdint.h>
#include "coblivious_layout.hpp"

typedef struct kd_s{
    unsigned int ID;
    unsigned char count;
} kd_t;


int test_oblivious(){
    size_t m = 1, beg = 32,
        end = 65;
    for(size_t i = beg; i  < end; i += m)
    {
        size_t n = i;

        std::vector<unsigned int> idt;
        std::vector<char> ct;
        std::vector<kd_t> tst;

        tst.resize(n);
        idt.resize(n);
        ct.resize(n);

        for(size_t j = 0u; j < n; j++) {
            tst[j].ID = j+1;
            tst[j].count = 50 + j;
        }


        ECDataCOLayout< std::vector<kd_t>::iterator, unsigned int, char >
            colayout(tst.begin(), n, 1);
        write_stvector(colayout.treeSteps, std::cout);
        write_stvector(colayout.stepPrefixes, std::cout);
        write_stvector(colayout.stepIds, std::cout);
        write_stvector(colayout.mIds, std::cout);
        write_stvector(colayout.mCounts, std::cout);
        char xCount;
        std::cout << colayout.getCount(n + 2, xCount) << " " << n + 2 << " ";
        std::cout << (int)xCount << ":";
        std::cout << colayout.getCount(n + 1, xCount) << " " << n + 1 << " ";
        std::cout << (int)xCount << ":";
        std::cout << colayout.getCount(0, xCount) << " " << 0 << " ";
        std::cout << (int)xCount << ":";
        std::cout << colayout.getCount(1, xCount) << " " << 1 << " ";
        std::cout << (int)xCount << ":";
        std::cout << colayout.getCount(31, xCount) << " " << 31 << " ";
        std::cout << (int)xCount << ":";
        std::cout << colayout.getCount(32, xCount) << " " << 32 << " ";
        std::cout << (int)xCount << ":";
        std::cout << colayout.getCount(33, xCount) << " " << 33 << " ";
        std::cout << (int)xCount << ":";
        std::cout << colayout.getCount(46, xCount) << " "  << 46 << " ";
        std::cout << (int)xCount << ":";
        std::cout << colayout.getCount(47, xCount) << " "  << 47 << " ";
        std::cout << (int)xCount << ":";
        std::cout << colayout.getCount(48, xCount) << " " << 48 << " ";
        std::cout << (int)xCount << ":";
        std::cout << std::endl;
        // std::cout << colayout.find(32)
        //           << std::endl;
    }
    return 0;

}
