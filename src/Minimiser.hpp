#ifndef SK7_MINIMISER_HPP
#define SK7_MINIMISER_HPP

#include "Kmer.hpp"

class Minimiser {
protected:
    uint64_t value;
    ushort size;
    uint64_t (*hashFunction) (Kmer, ushort);
public:
    Minimiser(uint64_t (*hashFunction) (Kmer, ushort), const Kmer& kmer, ushort size);
    uint64_t getValue() const;
    std::string toString() const;
};


#endif //SK7_MINIMISER_HPP
