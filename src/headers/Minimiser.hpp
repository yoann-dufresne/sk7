#ifndef SK7_MINIMISER_HPP
#define SK7_MINIMISER_HPP

#include "Kmer.hpp"

typedef struct hashPos {
    uint64_t hashValue;
    short pos;
} hashPos;

class Minimiser {
protected:
    uint64_t value;
    ushort length;
    short pos;
    hashPos (*hashFunction) (Kmer, ushort);

public:
    Minimiser(hashPos (*hashFunction) (Kmer, ushort), ushort length);
    Minimiser(hashPos (*hashFunction) (Kmer, ushort), ushort length, Kmer kmer);
    void init(Kmer kmer);
    void fromNewEnd(Kmer kmer);
    uint64_t getValue() const;
    short getPos() const;
    std::string toString() const;
};


#endif //SK7_MINIMISER_HPP
