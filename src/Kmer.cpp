#include "headers/Kmer.hpp"
#include <string>


/**
 * Constructor from value
 * @param value the value of the wanted kmer
 * @param length the length (in letter 'A', 'C', 'G', 'T') of the kmer
 */
Kmer::Kmer(uint64_t value, ushort length) {
    this->value = value;
    this->length = length;
}

/**
 * getter for the value
 * @return the value of the kmer
 */
uint64_t Kmer::getValue() {
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
                res = 'G' + res;
                break;
            default:
                res = 'T' + res;
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
ushort Kmer::getLength() {
    return this->length;
}

/**
 * get a subkmer of the original kmer with (0 <= start <= end < k)
 * @param start index of the wanted start of the sub kmer
 * @param end index of the end of the wanted kmer.
 * @return the Kmer[start:end+1]
 */
Kmer Kmer::getSubKmer(int start, int end) {
    uint64_t mask = (1 << ((this->length - start) * 2)) - 1;
    int decalage = ((this->length - end - 1) * 2);
    mask = mask >> decalage;
    mask = mask << decalage;
    return Kmer((this->value & mask) >> decalage, end - start + 1);
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
    return Kmer(((this->value >> ((this->length - pos) * 2)) << (suffixLen * 2)) + suffix, this->length - fragLength);
}
