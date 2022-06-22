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
    ushort length;
    uint64_t getValue();
    Kmer(uint64_t value, ushort length);
    std::string toString();
    ushort getLength();
    Kmer getSubKmer(int start, int end);
    Kmer removePart(int pos, int fragLength);
};


#endif //SK7_KMER_HPP
