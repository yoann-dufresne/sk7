#ifndef SK7_KMER_H
#define SK7_KMER_H

#include <cinttypes>

class Kmer {

protected:
    uint64_t value;
    ushort size;

public:
    uint64_t getValue();
    Kmer(uint64_t value, ushort size);
    ~Kmer() {};
};


#endif //SK7_KMER_H
