#ifndef SK7_KMER_H
#define SK7_KMER_H

#include <cinttypes>
#include <string>

class Kmer {

protected:
    uint64_t value;
    ushort size;

public:
    uint64_t getValue();
    Kmer(uint64_t value, ushort size);
    ~Kmer() = default;;
    std::string toString();
    ushort getSize();
};


#endif //SK7_KMER_H
