#include <vector>
#include <iostream>
#include <iterator>
#include <cassert>
#include <cstdint>
#include "implicit_heap_search_c.h"

typedef std::vector<uint64_t>::size_type size_type;
typedef std::vector<size_type> size_type_vector;

template<typename T>
void write_stvector(std::vector<T> sts, std::ostream& os) {
    std::ostream_iterator<T> outs(os, " ");
    std::copy(sts.cbegin(), sts.cend(), outs);
    os << std::endl;
}


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

class CacheAwareLayoutHelper{
   private:
    size_type_vector completeLevelSize;
    size_type_vector completeTreeSize;
    size_type d_cache = 0;
    size_type k_cache = 0;
    size_type m_cache = 0;
    size_type_vector stsize;

    void init_sizes(const size_type &m,
                    const size_type &l){
        completeLevelSize.resize(l);
        completeTreeSize.resize(l);
        completeTreeSize[0] =  completeLevelSize[0] = m;
        for(auto i = 1u; i < l; i++){
            completeLevelSize[i] = completeLevelSize[i - 1] * (m + 1);
            completeTreeSize[i] = completeTreeSize[i - 1] + completeLevelSize[i];
        }
#ifdef NDEBUG
    write_stvector(completeLevelSize, std::cout);
    write_stvector(completeTreeSize, std::cout);
#endif
    }


  public:
    CacheAwareLayoutHelper(const size_type& n,
                           const size_type& l,
                           const size_type& m):d_cache(n),k_cache(l),m_cache(m){
        init_sizes(m,l);
        subtree_size(n, l, m);
    }

    size_type_vector& subtree_size(const size_type& d, const size_type& k,
                                   const size_type& m){
        size_type_vector ls;
        stsize.resize(m + 1);
        ls.resize(m + 1);
        // No. of elements at last level.
        auto dk = d - completeTreeSize[k - 2];
        //
        auto q = dk / completeLevelSize[k - 2];
        auto r = dk % completeLevelSize[k - 2];
#ifdef NDEBUG
        std::cout << "dk : " << dk
                  << " q : " << q
                  << " r : " << r
                  << std::endl;
#endif
        //TODO: cache
        for(auto j = 0u; j <= m; j++)
            if(j < q)
                ls[j] = completeLevelSize[k - 2];
            else if(j > q)
                ls[j] = 0;
            else
                ls[j] = r;

#ifdef NDEBUG
        std::cout << "LS : ";
        write_stvector(ls, std::cout);
#endif
        //TODO: cache st size
        for(auto j = 0u; j <= m; j++)
            if(k > 2)
                stsize[j] = ls[j] + completeTreeSize[k - 3];
            else
                stsize[j] = ls[j];
#ifdef NDEBUG
        std::cout << "ST : ";
        write_stvector(stsize, std::cout);
#endif
        return stsize;
    }

    size_type_vector& subtree_size_cache(const size_type& d,
            const size_type& k,
            const size_type& m){
        if(d_cache == d && m_cache == m && k_cache == k)
            return stsize;
        // initialization
        d_cache = d;
        k_cache = k;
        m_cache = m;
        return subtree_size(d,k,m);
    }
};


template <typename T>
void caware_layout(std::vector<T>& srtx,
                   const size_type& m,
                   std::vector<T>& ctree){
    auto n = srtx.size();
    assert(m > 2);
    assert(n > 0);
    assert(n % m == 0);
    // total no. of levels
    size_type l = int_log(m + 1, n + 1);
#ifdef NDEBUG
    std::cout << "Levels : " << l << std::endl;
#endif

    CacheAwareLayoutHelper chelper(n, l, m);
    size_type i = 0, l_ptr = 0, c_ptr = 0;

    // Root node indicies and level
    ctree[0] = 0;
    ctree[1] = l;
    ctree[m - 1] = n - 1;
    while(i < n){
        // [x,y] is the range of indicies covered by this sub-tree
        auto x = ctree[c_ptr];
        auto y = ctree[c_ptr + m - 1];
        // k no. levels of current sub-tree including root
        size_type k = ctree[c_ptr + 1];
        size_type d = y - x + 1;
#ifdef NDEBUG
        std::cout << "x : " << x
                  << " y : " << y
                  << " k : " << k
                  << std::endl;
#endif
        if(d == m){
            // no sub-trees
            for(auto j = 0u; j < m; j++)
                ctree[c_ptr + j] = srtx[x++];
        } else {
            const size_type_vector& stsize =
                chelper.subtree_size_cache(d,k,m);

            y = x;
            // Update entries for current node
            for(auto j = 0u; j < m; j++) {
                y = y + stsize[j];
                ctree[c_ptr + j] = srtx[y];
                y++;
            }
            // Insert indices for sub-trees
            for(auto j = 0u; j <= m; j++) {
                if(stsize[j] > 0) {
                    y = x + stsize[j];
                    l_ptr = l_ptr + m;
                    ctree[l_ptr] = x;
                    ctree[l_ptr + 1] = k - 1;
                    ctree[l_ptr + m - 1] = y - 1;
                    x = y + 1;
                }
            }
        }
        // move current forward
        i += m;
        c_ptr += m;
    }
}

int main(int argc, char* argv[]){

    size_type m = 3, beg = m * 8,
        end = 1 + (m * 40);
    for(size_type i = beg; i  < end; i += m) {
        size_type_vector stsize;
        size_type n = i;

        std::vector<int> tst,ct;
        tst.resize(n);
        ct.resize(n);

        for(auto i = 0u; i < n; i++)
            tst[i] = i+1;

        caware_layout(tst, m, ct);
        write_stvector(ct, std::cout);

        std::cout << "Search : ";
        for(auto j = -1; j <= (int) (i+2);j++) {
            std::cout << pure_c::implicit_heap_search(ct.begin(), ct.end(), m + 1, j) << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}
