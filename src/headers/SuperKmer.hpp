#ifndef SK7_SUPERKMER_HPP
#define SK7_SUPERKMER_HPP

#include <cinttypes>
#include <vector>
#include <set>

#include "Kmer.hpp"


#define TYPE uint8_t
#define SIZE 8 // For uint8_t

class SuperKmer {

public:
    /// Attributes
    std::vector<TYPE> tab;
    enum logic {SUPERIOR, INFERIOR, EQUAL, INCOMPARABLE};

    /// Constructor
    SuperKmer(std::vector<TYPE> tab);
    SuperKmer();

    /// Access and modification
    uint64_t accessBits(int start, int end);
    void setBits(const int &start, const int &length, const uint64_t &value);

    /// Getter
    uint64_t getValue();
    int getPrefixLen();
    int getSuffixLen();

    /// Storage
    std::vector<SuperKmer> cut(const int &commonPartStart, const int &commonPartEnd, const int &fixBitSize); // To redo to cut specific part
    static std::vector<logic> compareSK(SuperKmer& superKmer1, SuperKmer& superKmer2);
    uint64_t buildSKMask(const int &prefixLen, const int &suffixLen);

    /// Operator
    bool operator==(const SuperKmer &toCompare) const;
    SuperKmer operator&(SuperKmer toIntersect); // Intersection
    std::vector<SuperKmer> operator^(SuperKmer &toXor); // Xor

    /// Intern methods
    std::vector<SuperKmer> split();
    uint64_t nonInterleavedKmerValue();

    /// Misc.
    void print();

};

#endif //SK7_SUPERKMER_HPP
