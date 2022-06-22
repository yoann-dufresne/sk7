#ifndef SK7_SUPERKMER_HPP
#define SK7_SUPERKMER_HPP

#include <cinttypes>
#include <vector>

#include "Kmer.hpp"


#define TYPE uint8_t
#define SIZE 8 // For uint8_t

class SuperKmer {

public:
    /// Attributes
    std::vector<TYPE> tab;

    /// Constructor
    SuperKmer(std::vector<TYPE> tab);
    SuperKmer();

    /// Access and modification
    uint64_t accessBits(int start, int end);
    void setBits(const int &start, const int &length, const uint64_t &value);

    /// Getter
    uint64_t getValue(const int &fixBitSize);
    int getPrefixLen(const int &fixBitSize);
    int getSuffixLen(const int &fixBitSize);

    /// Storage
    std::vector<SuperKmer> cut(const int &commonPartStart, const int &commonPartEnd, const int &fixBitSize);

    /// Operator
    bool operator==(SuperKmer toCompare);
    SuperKmer operator&(SuperKmer toIntersect); // Intersection

    /// Misc.
    void print();

};

#endif //SK7_SUPERKMER_HPP
