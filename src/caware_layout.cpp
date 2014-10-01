#include <vector>
#include <iostream>
#include <iterator>
#include <cassert>
#include <stdint.h>
#include "implicit_heap_search_c.h"
#include "caware_layout.hpp"

size_type int_log(const size_type& b,
                  const size_type& x){
    double st = b;
    size_type vl = 1;
    while(st < x){
        st = st * b;
        vl++;
    }
    return vl;
}


typedef struct kd_s{
    int ID;
    unsigned char count;
} kd_t;

int test_caware(){

    size_type m = 3, beg = m * 8,
        end = 1 + (m * 40);
    for(size_type i = beg; i  < end; i += m) {
        size_type_vector stsize;
        size_type n = i;

        std::vector<int> idt;
        std::vector<char> ct;
        std::vector<kd_t> tst;

        tst.resize(n);
        idt.resize(n);
        ct.resize(n);

        for(size_type j = 0u; j < n; j++) {
            tst[j].ID = j+1;
            tst[j].count = n - j;
        }

        caware_layout(tst.begin(), tst.end(), m, idt, ct);
        write_stvector(idt, std::cout);

        std::cout << "Search : ";
        for(int j = -1; j <= (int) (i+2);j++) {
            std::vector<int>::iterator where = pure_c::implicit_heap_search(idt.begin(), idt.end(),
                                                      m + 1, j);
            std::cout << (where != idt.end()) << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}
