#ifndef _IO_UTIL_H_
#define _IO_UTIL_H_
#include <vector>
#include <iterator>
#include <iostream>

template<typename T>
void write_stvector(std::vector<T> sts, std::ostream& os) {
    std::ostream_iterator<T> outs(os, " ");
    std::copy(sts.begin(), sts.end(), outs);
    os << std::endl;
}


#endif /* _IO_UTIL_H_ */
