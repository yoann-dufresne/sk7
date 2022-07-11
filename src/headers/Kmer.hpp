#ifndef SK7_KMER_HPP
#define SK7_KMER_HPP

#include <cinttypes>
#include <string>
#include <cmath>

namespace sk7 {
    extern int k;
    extern int m;
    extern int fixBitSize;
    void initLib(int _k, int _m);
}

class Kmer {

protected:
    uint64_t value;

public:
    /// Attributes
    ushort length;

    /// Constructor
    Kmer();
    Kmer(uint64_t value);
    Kmer(uint64_t value, ushort length);

    /// Getter
    uint64_t getValue();
    ushort getLength();

    /// Utils
    Kmer getSubKmer(int start, int end);
    Kmer removePart(int pos, int fragLength);
    Kmer reverseComplement();

    /// Operator
    bool operator<(const Kmer &toCompare) const;
    bool operator==(const Kmer &toCompare) const;

    /// Misc.
    std::string toString();

};


#endif //SK7_KMER_HPP
