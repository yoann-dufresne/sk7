#ifndef SK7_KMER_HPP
#define SK7_KMER_HPP

#include <cinttypes>
#include <string>
#include <cmath>

#include "utils.hpp"

/// namespace for share use of variables or initialisation
namespace sk7 {

typedef unsigned __int128 uint128_t;

class Kmer {

    protected:
        uint128_t value;

    public:
        /// Attributes
        ushort length;

        /// Constructor
        Kmer();
        Kmer(uint128_t value);
        Kmer(uint128_t value, ushort length);

        /// Getter
        uint128_t getValue() const;
        ushort getLength() const;

        /// Utils
        Kmer getSubKmer(int start, int end);
        Kmer removePart(int pos, int fragLength);
        Kmer reverseComplement();

        /// Operator / Comparator
        bool operator<(const Kmer &toCompare) const;
        bool operator==(const Kmer &toCompare) const;
        bool operator!=(const Kmer &toCompare) const;
        bool operator>(const Kmer &toCompare) const;
        bool operator>=(const Kmer &toCompare) const;
        static bool fullComparison(const Kmer &, const Kmer &);

        /// Misc.
        std::string toString();

    };


    extern int k;
    extern int m;
    extern int fixBitSize;

    extern bool (*infKmer) (const Kmer &, const Kmer &);
    extern bool (*equalKmer) (const Kmer &, const Kmer &);

    void initLib(int _k, int _m, bool (*_infKmer) (const Kmer &, const Kmer &) = &infId, bool (*_equalKmer) (const Kmer &, const Kmer &) = &equalId);
}

#endif //SK7_KMER_HPP
