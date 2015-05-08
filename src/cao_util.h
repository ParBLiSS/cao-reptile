#ifndef _CAO_UTIL_H_
#define _CAO_UTIL_H_

#include <vector>
#include <iterator>
#include <iostream>

typedef std::vector<uint64_t>::size_type size_type;
typedef std::vector<size_type> size_type_vector;

template<typename T>
void write_stvector(std::vector<T> sts, std::ostream& os) {
    std::ostream_iterator<T> outs(os, " ");
    std::copy(sts.begin(), sts.end(), outs);
    os << std::endl;
}


#endif /* _CAO_UTIL_H_ */
