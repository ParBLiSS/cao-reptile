#ifndef _CAWARE_LAYOUT_H_
#define _CAWARE_LAYOUT_H_
#include <vector>
#include <iostream>
#include <iterator>
#include <cassert>
#include <cstdint>
#include <fstream>

#include "aligned_allocator.hpp"
#include "io_util.h"
#include "cao_util.h"
#include "implicit_heap_search_c.h"

size_type int_log(const size_type& b, const size_type& x);
int test_caware();

class CacheAwareLayoutHelper{
   private:
    size_type_vector completeLevelSize;
    size_type_vector completeTreeSize;
    size_type d_cache;
    size_type k_cache;
    size_type m_cache;
    size_type_vector stsize;

    void init_sizes(const size_type &m,
                    const size_type &l){
        completeLevelSize.resize(l);
        completeTreeSize.resize(l);
        completeTreeSize[0] =  completeLevelSize[0] = m;
        for(size_type i = 1u; i < l; i++){
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
        size_type dk = d - completeTreeSize[k - 2];
        //
        size_type q = dk / completeLevelSize[k - 2];
        size_type r = dk % completeLevelSize[k - 2];
#ifdef NDEBUG
        std::cout << "dk : " << dk
                  << " q : " << q
                  << " r : " << r
                  << std::endl;
#endif
        //TODO: cache
        for(size_type j = 0u; j <= m; j++)
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
        for(size_type j = 0u; j <= m; j++)
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


template <typename RandomIterator,
          typename IdIterator,
          typename CountIterator>
void caware_layout(RandomIterator srt_begin,
                   RandomIterator srt_end,
                   const size_type& m,
                   IdIterator id_tree,
                   CountIterator count_tree){

    size_type n = std::distance(srt_begin, srt_end);
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
    id_tree[0] = 0;
    id_tree[1] = l;
    id_tree[m - 1] = n - 1;
    while(i < n){
        // [x,y] is the range of indicies covered by this sub-tree
        size_type x = id_tree[c_ptr];
        size_type y = id_tree[c_ptr + m - 1];
        // k no. levels of current sub-tree including root
        size_type k = id_tree[c_ptr + 1];
        size_type d = y - x + 1;
#ifdef NDEBUG
        std::cout << "x : " << x
                  << " y : " << y
                  << " k : " << k
                  << std::endl;
#endif
        if(d == m){
            // no sub-trees
            for(size_type j = 0u; j < m; j++) {
                id_tree[c_ptr + j] = (srt_begin + x)->ID;
                count_tree[c_ptr + j] = (srt_begin + x)->count;
                x++;
            }
        } else {
            const size_type_vector& stsize =
                chelper.subtree_size_cache(d,k,m);

            y = x;
            // Update entries for current node
            for(size_type j = 0u; j < m; j++) {
                y = y + stsize[j];
                id_tree[c_ptr + j] = (srt_begin + y)->ID;
                count_tree[c_ptr + j] = (srt_begin + y)->count;
                //id_tree[c_ptr + j] = srtx[y];
                y++;
            }
            // Insert indices for sub-trees
            for(size_type j = 0u; j <= m; j++) {
                if(stsize[j] > 0) {
                    y = x + stsize[j];
                    l_ptr = l_ptr + m;
                    id_tree[l_ptr] = x;
                    id_tree[l_ptr + 1] = k - 1;
                    id_tree[l_ptr + m - 1] = y - 1;
                    x = y + 1;
                }
            }
        }
        // move current forward
        i += m;
        c_ptr += m;
    }
}

template <typename RandomIterator,
          typename IdType,
          typename CountType>
class ECDataCALayout{

    int nbytes(IdType val){
        int n = 0;
        if(val == 0)
            return 1;
        while (val != 0) {
            val >>= 8;
            n ++;
        }
        return n;
    }

public:
    typedef typename std::vector<IdType, aligned_allocator<IdType, 64> >::iterator id_iterator;
        typedef typename std::vector<IdType, aligned_allocator<IdType, 64> >::const_iterator id_const_iterator;
    std::vector< IdType, aligned_allocator<IdType, 64> > mIds;
    std::vector<CountType, aligned_allocator<CountType, 64> > mCounts;
    int mDegree;
    ECDataCALayout(){};
    void init(RandomIterator intr,
              size_type nTotal, int kSize){
        assert(nTotal > 0);
        assert(nTotal % kSize == 0);
        mIds.resize(nTotal);
        mCounts.resize(nTotal);
        mDegree = kSize + 1;
        caware_layout(intr, intr + nTotal,
                      kSize, mIds.begin(), mCounts.begin());
    };

    ECDataCALayout(RandomIterator intr,
                 size_type nTotal,
                 int kSize){
        init_layout(intr, nTotal, kSize);
    };

    bool find(const IdType &rID) const{
        return (pure_c::implicit_heap_search(mIds.begin(),
                                             mIds.end(),
                                             mDegree,
                                             rID)
                !=  mIds.end());
    };

    int getCount(const IdType& rID, CountType& count) const{
        id_const_iterator fpos = pure_c::implicit_heap_search(mIds.begin(),
                                                              mIds.end(),
                                                              mDegree,
                                                              rID);
        if(fpos != mIds.end()) {
            int dist = fpos - mIds.begin();
            count = mCounts[dist];
            //std::cout << "count " << (int) count << std::endl;
            return dist;
        } else {
            return -1;
        }
    };

    void compression_stats(double& avg_bytes,
                           double& min_bytes,
                           double& max_bytes){
        unsigned m = mDegree - 1;
        IdType e0, d0;
        std::vector<IdType> cline;
        cline.resize(m);
        avg_bytes = max_bytes = 0.0;
        min_bytes = sizeof(IdType) * m;
        for(unsigned i = 0; i < mIds.size(); i += m){
            double tbytes = 0;
            e0 = mIds[i];
            for(unsigned j = 1;  j < m; j++)
                cline[j] = mIds[i + j] - mIds[i];
            d0 = cline[1];
            for(unsigned j = 2;  j < m; j++)
                cline[j] = cline[j] - cline[1];
            tbytes += nbytes(e0);
            tbytes += nbytes(d0);
            for(unsigned j = 2;  j < m; j++)
                tbytes += nbytes(cline[j]);
            if(tbytes > max_bytes)
                max_bytes = tbytes;
            if(tbytes < min_bytes)
                min_bytes = tbytes;
            avg_bytes += tbytes;
        }
        avg_bytes = (avg_bytes * m)/mIds.size();
    }

    void serialize(std::ofstream& ofs) const{
        ofs.write((char*)&mIds[0], sizeof(IdType)*mIds.size());
        ofs.write((char*)&mCounts[0], sizeof(CountType)*mCounts.size());
    }
};


#endif /* _CAWARE_LAYOUT_H_ */
