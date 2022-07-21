#ifndef SK7_UTILS_HPP
#define SK7_UTILS_HPP

#include <iostream>

class Kmer;

/// Kmer manipulation
uint64_t interleavedOrder(Kmer &kmer, int minimiserPos);
uint64_t reorderValue(uint64_t value, const int &prefixLen, const int &suffixLen);

/// Example for Kmer comparison
bool infId(const Kmer &kmer1, const Kmer &kmer2);
bool equalId(const Kmer &kmer1, const Kmer &kmer2);

#include "Kmer.hpp"

#endif //SK7_UTILS_HPP
