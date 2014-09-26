#ifndef _IO_UTIL_H_
#define _IO_UTIL_H_


template<typename T>
void write_stvector(std::vector<T> sts, std::ostream& os) {
    std::ostream_iterator<T> outs(os, " ");
    std::copy(sts.cbegin(), sts.cend(), outs);
    os << std::endl;
}


#endif /* _IO_UTIL_H_ */
