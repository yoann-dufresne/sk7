#ifndef SK7_SUPERKMER_HPP
#define SK7_SUPERKMER_HPP

#include <cinttypes>
#include <vector>
#include <bitset> // for tests


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
    uint64_t accessBits(int start, int end) const;
    void setBits(const int &start, const int &length, const uint64_t &value);

    /// Getter
    uint64_t getValue() const;
    int getPrefixLen() const;
    int getSuffixLen() const;

    /// Storage
    std::vector<SuperKmer> cut(const int &commonPartStart, const int &commonPartEnd, const int &fixBitSize); // To redo to cut specific part
    static std::vector<logic> compareSK(SuperKmer& superKmer1, SuperKmer& superKmer2);
    uint64_t buildSKMask(const int &prefixLen, const int &suffixLen) const;
    Kmer readKmer(int kmerPrefixLen) const;
    SuperKmer extract(int wantedPrefixLen) const;

    /// Operator
    bool operator==(const SuperKmer &toCompare) const;
    SuperKmer operator&(SuperKmer toIntersect); // Intersection
    SuperKmer operator|(SuperKmer toUnite); // Union

    /// Intern methods
    std::vector<SuperKmer> split() const;
    uint64_t nonInterleavedKmerValue() const;

    /// Misc.
    void print() const;

};

#endif //SK7_SUPERKMER_HPP
