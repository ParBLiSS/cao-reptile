#include <vector>
#include <iostream>
#include <iterator>
#include <cassert>
#include <stdint.h>
#include "coblivious_layout.hpp"
#include "build_hp_tree.hpp"
#include "hp_precomputed_table.h"
#include "implicit_hp_search_pure_c.hpp"
#include "io_util.h"
#include "log2.h"

typedef struct kd_s{
    int ID;
    unsigned char count;
} kd_t;

int test_oblivious(){
    size_t m = 4, beg = m * 16,
        end = 1 + (m * 16), i = 64;
    //for(size_t i = beg; i  < end; i += m)
    {
        size_t n = i;

        std::vector<int> idt;
        std::vector<char> ct;
        std::vector<kd_t> tst;

        tst.resize(n);
        idt.resize(n);
        ct.resize(n);

        for(size_t i = 0u; i < n; i++) {
            tst[i].ID = i+1;
            tst[i].count = n - i;
        }

        coblivious_layout(tst.begin(), tst.end(), idt.begin(),
                          idt.end(),ct.begin(), ct.end(), 6, 1);
        write_stvector(idt, std::cout);
        // idt.resize(1);
        // ct.resize(1);
        // idt.resize(n);
        // ct.resize(n);
        // coblivious_layout(tst.begin(), tst.end(), idt.begin(),
        //                   idt.end(),ct.begin(), ct.end(), 3, 1);
        // write_stvector(idt, std::cout);
        // idt.resize(1);
        // ct.resize(1);
        // idt.resize(n);
        // ct.resize(n);
        // coblivious_layout(tst.begin(), tst.end(), idt.begin(),
        //                   idt.end(),ct.begin(), ct.end(), 2, 1);
        // write_stvector(idt, std::cout);
        precomputed_table<int> pct(6);

        std::cout << pure_c::implicit_hp_search(idt.begin(), idt.end(),
                                                pct.Pos.begin(),
                                                pct.T.begin(),
                                                pct.B.begin(),
                                                pct.D.begin(), 0)
                  << " "
                  << pure_c::implicit_hp_search(idt.begin(), idt.end(),
                                                pct.Pos.begin(),
                                                pct.T.begin(),
                                                pct.B.begin(),
                                                pct.D.begin(), 61)
                  << std::endl;
        // std::cout << "Search : ";
        // for(int j = -1; j <= (int) (i+2);j++) {
        //     auto where = pure_c::implicit_heap_search(idt.begin(), idt.end(),
        //                                               m + 1, j);
        //     std::cout << (where != idt.end()) << " ";
        // }
        // std::cout << std::endl;
    }
    return 0;

}
