#ifndef SK7_UTILS_HPP
#define SK7_UTILS_HPP

#include "Kmer.hpp"

#include <iostream>

uint64_t interleavedOrder(Kmer &kmer, int minimiserPos);
uint64_t reorderValue(uint64_t value, const int &prefixLen, const int &suffixLen);

#endif //SK7_UTILS_HPP
