#ifndef ALIGNED_ALLOCATOR_H
#define ALIGNED_ALLOCATOR_H
#include <stdlib.h>
#include <stdexcept>
#include <cstddef>
// Mainly from
// https://gist.githubusercontent.com/donny-dont/1471329/raw/8f063f5f4326d301b14fe6b781495be54ca48941/aligned_allocator.cpp
// and
// http://stackoverflow.com/questions/17378444/stdalign-and-stdaligned-storage-for-aligned-allocation-of-memory-blocks

// void* allocate_aligned_memory(size_t align, size_t size){
//     assert(align >= sizeof(void*));

//     if (size == 0)
//         return (void*)0;

//     void* ptr = 0;
//     int rc = posix_memalign(&ptr, align, size);

//     return ptr;
// }

// void deallocate_aligned_memory(void *ptr){
//     return free(ptr);
// }

template <typename T, std::size_t Alignment>
class aligned_allocator {
public:
    typedef T * pointer;
    typedef const T * const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T value_type;
    typedef std::size_t size_type;
    typedef ptrdiff_t difference_type;

    T * address(T& r) const {
        return &r;
    }

    const T * address(const T& s) const{
        return &s;
    }

    std::size_t max_size() const{
        return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
    }

    template <typename U>
    struct rebind {
            typedef aligned_allocator<U, Alignment> other;
    } ;

    bool operator!=(const aligned_allocator& other) const {
        return !(*this == other);
    }

    void construct(T * const p, const T& t) const {
        void * const pv = static_cast<void *>(p);
        new (pv) T(t);
    }

    void destroy(T * const p) const {
        p->~T();
    }

    bool operator==(const aligned_allocator&) const{
        return true;
    }

    aligned_allocator() { }
    aligned_allocator(const aligned_allocator&) { }
    template <typename U> aligned_allocator(const aligned_allocator<U, Alignment>&) { }
    ~aligned_allocator() { }

    T * allocate(const std::size_t n) const{
        if (n == 0)
            return NULL;
        if (n > max_size())
            throw std::length_error(
                    "aligned_allocator<T>::allocate() - Integer overflow.");
        void *pv;
        int rc = posix_memalign(&pv, Alignment, n * sizeof(T));
        if (rc != 0 || pv == NULL)
            throw std::bad_alloc();
        return static_cast<T *>(pv);
    }

    void deallocate(T * const p, const std::size_t) const{
        free(p);
    }

    template <typename U>
    T* allocate(const std::size_t n, const U * /* const hint */) const{
        return allocate(n);
    }

private:
    aligned_allocator& operator=(const aligned_allocator&);
};


#endif /* ALIGNED_ALLOCATOR_H */
