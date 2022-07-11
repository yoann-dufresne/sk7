#ifndef SK7_MINIMISER_HPP
#define SK7_MINIMISER_HPP

#include "Kmer.hpp"
#include "utils.hpp"

typedef struct hashPos {
    uint64_t hashValue;
    short pos;
} hashPos;

class Minimiser {
protected:
    /// Attributes
    uint64_t value;
    ushort length;
    short pos;
    hashPos (*hashFunction) (Kmer, ushort);

public:
    /// Constructor and related
    Minimiser(hashPos (*hashFunction) (Kmer, ushort), ushort length);
    Minimiser(hashPos (*hashFunction) (Kmer, ushort), ushort length, Kmer kmer);
    Minimiser(hashPos (*hashFunction) (Kmer, ushort), Kmer &kmer);
    void init(Kmer kmer);
    void fromNewEnd(Kmer kmer);

    /// Getter
    uint64_t getValue() const;
    short getPos() const;

    /// Misc.
    std::string toString() const;
};


#endif //SK7_MINIMISER_HPP
