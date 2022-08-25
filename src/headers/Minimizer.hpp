#ifndef SK7_MINIMIZER_HPP
#define SK7_MINIMIZER_HPP

#include "Kmer.hpp"
#include "utils.hpp"

namespace sk7 {

typedef struct hashPos {
    uint64_t hashValue;
    short pos;
} hashPos;

class Minimizer {
protected:
    /// Attributes
    uint64_t value;
    ushort length;
    short pos;

    hashPos (*hashFunction)(Kmer, ushort);

public:
    /// Constructor
    Minimizer(hashPos (*hashFunction)(Kmer, ushort), ushort length, Kmer kmer);
    Minimizer(hashPos (*hashFunction)(Kmer, ushort), Kmer &kmer);
    /// Getter
    uint64_t getValue() const;
    short getPos() const;

    /// Misc.
    std::string toString() const;
};

}
#endif //SK7_MINIMIZER_HPP
