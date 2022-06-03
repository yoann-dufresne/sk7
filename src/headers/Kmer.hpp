#ifndef SK7_KMER_H
#define SK7_KMER_H

#include <cinttypes>
#include <string>

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


#endif //SK7_KMER_H
