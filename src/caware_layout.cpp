#include <vector>
#include <iostream>
#include <iterator>
#include <cassert>
#include <stdint.h>
#include "implicit_heap_search_c.h"
#include "caware_layout.hpp"



typedef struct kd_s{
    int ID;
    unsigned char count;
} kd_t;

int test_caware(){
    typedef std::vector<kd_t>::iterator kd_iterator;
    size_type m = 3, beg = m * 8,
        end = 1 + (m * 40), nt = 2;
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
            tst[j].count = 34+ (j % 94);
        }

        caware_layout(tst.begin(), tst.end(), m, idt.begin(), ct.begin());
        write_stvector(idt, std::cout);

        std::cout << "Search : ";
        for(int j = -1; j <= (int) (i+2);j++) {
            std::vector<int>::iterator where = pure_c::implicit_heap_search(idt.begin(), idt.end(),
                                                      m + 1, j);
            std::cout << (where != idt.end()) << " ";
        }
        std::cout << std::endl;
        //ECDataCALayout<kd_iterator, int, char> eca(tst.begin(), tst.size(), m);

        while(tst.size() % (m * nt) > 0)
            tst.push_back(tst.back());
        // std::cout << "T:" << tst.size() << std::endl;
        ParECDataCALayout<kd_iterator, int, char> pca(tst.begin(), tst.size(),
                                                      m, nt);
        std::cout << "Search : ";
        for(int j = -1; j <= (int) (i+2);j++) {
            char cnt = 33;
            pca.getCount(j, cnt);
            std::cout << cnt << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}
