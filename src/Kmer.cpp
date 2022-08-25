#include "Kmer.hpp"

#include <string>

/**
 * A namespace used to store variables that we want to share across the classes
 */
namespace sk7 {
    int k;
    int m;
    int fixBitSize;

    bool (*infKmer)(const Kmer &, const Kmer &);

    bool (*equalKmer)(const Kmer &, const Kmer &);

    /**
     * Initialise the global parameters of the library
     * @param _k the size of the Kmers
     * @param _m the size of the Minimisers
     * @param _infKmer a function to order two Kmers (all functions are not guaranteed to work ex. those who count "A", compare specific nucleotide at specific indexes)
     * @param _equalKmer a function to verify the equality between two Kmers (same)
     */
    void initLib(int _k, int _m, bool (*_infKmer)(const Kmer &, const Kmer &),
                 bool (*_equalKmer)(const Kmer &, const Kmer &)) {
        k = _k;
        m = _m;
        fixBitSize = ceil(log2(k - m + 1));
        infKmer = _infKmer;
        equalKmer = _equalKmer;
    }


/**
 * Default constructor
 */
    Kmer::Kmer() {
        this->value = 0;
        this->length = sk7::k;
    }

/**
 * Constructor to build a Kmer from its value
 * @param value the value of the Kmer
 */
    Kmer::Kmer(uint128_t value) {
        this->value = value;
        this->length = sk7::k;
    }

/**
 * Constructor from value and length
 * @param value the value of the wanted kmer
 * @param length the length (in letter 'A', 'C', 'G', 'T') of the kmer
 */
    Kmer::Kmer(uint128_t value, ushort length) {
        this->value = value;
        this->length = length;
    }

/**
 * getter for the value
 * @return the value of the kmer
 */
    uint128_t Kmer::getValue() const {
        return this->value;
    }

/**
 * Give a textual representation of a kmer
 * @return a string on {'A', 'C', 'G', 'T'} that represents the kmer
 */
    std::string Kmer::toString() {
        uint64_t tmp = value;
        std::string res;
        res.reserve(this->length);
        ushort b = 4;
        while (tmp != 0) {
            uint64_t q = tmp / b;
            int r = tmp % b;
            switch (r) {
                case 0:
                    res = 'A' + res;
                    break;
                case 1:
                    res = 'C' + res;
                    break;
                case 2:
                    res = 'T' + res;
                    break;
                default:
                    res = 'G' + res;
                    break;
            }
            tmp = q;
        }
        while (res.length() < this->length) {
            res = "A" + res;
        }
        return res;
    }

/**
 * a getter for the Kmer length
 * @return the length of the kmer
 */
    ushort Kmer::getLength() const {
        return this->length;
    }

/**
 * get a subkmer of the original kmer with (0 <= start <= end < k)
 * @param start index of the wanted start of the sub kmer
 * @param end index of the end of the wanted kmer (included)
 * @return the Kmer from start to end (both included) of the calling kmer.
 */
    Kmer Kmer::getSubKmer(int start, int end) {
        uint64_t mask = (1 << ((this->length - start) * 2)) - 1;
        int shift = ((this->length - end - 1) * 2);
        // set the 0s of the mask
        mask = mask >> shift;
        mask = mask << shift;
        return Kmer((this->value & mask) >> shift, end - start + 1);
    }

/**
 * Remove a part from a Kmer (ex. "ACGTACGT".removePart(2, 3) = "ACCGT")
 * @param pos the starting position of the unwanted fragment
 * @param fragLength the fragLength of the unwanted fragment
 * @return a Kmer without the unwanted fragment
 */
    Kmer Kmer::removePart(int pos, int fragLength) {
        int suffixLen = this->length - pos - fragLength;
        uint64_t mask = (1 << ((suffixLen) * 2)) - 1;
        uint64_t suffix = this->value & mask;
        return Kmer(((this->value >> ((this->length - pos) * 2)) << (suffixLen * 2)) + suffix,
                    this->length - fragLength);
    }


/**
 * Build and return the reverse complement of this Kmer
 * @return a reverse Kmer
 */
    Kmer Kmer::reverseComplement() {
        uint64_t res = value;

        res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
        res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
        res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
        res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
        res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);

        res ^= 0xAAAAAAAAAAAAAAAA;
        res >>= (2 * (32 - sk7::k));

        return Kmer(res);
    }

/**
 * Compare the value of two Kmers
 * @param toCompare the Kmer to compare with
 * @return this < toCompare with a custom function if provided with InitLib
 */
    bool Kmer::operator<(const Kmer &toCompare) const {
        return sk7::infKmer(*this, toCompare);
    }

/**
 * Equality between two Kmers with a personalized function
 * @param ToCompare the Kmer to compare with
 * @return true if the Kmers are the same (for the given function), else false
 */
    bool Kmer::operator==(const Kmer &toCompare) const {
        return sk7::equalKmer(*this, toCompare);
    }

/**
 * Difference between two Kmers with a personalized function
 * @param ToCompare the Kmer to compare with
 * @return true if the Kmers are the same (for the given function), else false
 */
    bool Kmer::operator!=(const Kmer &toCompare) const {
        return not(*this == toCompare);
    }

/**
 * Compare the value of two Kmers
 * @param toCompare the Kmer to compare with
 * @return this > toCompare with a custom function if provided with InitLib
 */
    bool Kmer::operator>(const Kmer &toCompare) const {
        return not(*this == toCompare or *this < toCompare);
    }

/**
 * Compare the value of two Kmers
 * @param toCompare the Kmer to compare with
 * @return this >= toCompare with a custom function if provided with InitLib
 */
    bool Kmer::operator>=(const Kmer &toCompare) const {
        return not(*this < toCompare);
    }

/**
 * Compare the plain parameters of two Kmers
 * @param kmer1 the first Kmer
 * @param kmer2 the second Kmer
 * @return true if the Kmers are the same, else false
 */
    bool Kmer::fullComparison(const Kmer &kmer1, const Kmer &kmer2) {
        return kmer1.length == kmer2.length && kmer1.value == kmer2.value;
    }


}