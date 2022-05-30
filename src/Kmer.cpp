#include <string>

#include "headers/Kmer.hpp"


/**
 * Constructor from value
 * @param value the value of the wanted kmer
 * @param size the size (in letter 'A', 'C', 'G', 'T') of the kmer
 */
Kmer::Kmer(uint64_t value, ushort size) {
    this->value = value;
    this->size = size;
}

/**
 * getter for the value
 * @return the value of the kmer
 */
uint64_t Kmer::getValue() {
    return this->value;
}

/**
 * Give a textual reprensentation of a kmer
 * @return a string on {'A', 'C', 'G', 'T'} that represents the kmer
 */
std::string Kmer::toString() {
    uint64_t tmp = value;
    std::string res;
    ushort b = 4;
    while (tmp != 0) {
        uint64_t q = tmp / b;
        int r = tmp % b;
        switch (r) {
            case 0:
                res = "A" + res;
                break;
            case 1:
                res = "C" + res;
                break;
            case 2:
                res = "G" + res;
                break;
            default:
                res = "T" + res;
        }
        tmp = q;
    }
    while (res.length() < this->size) {
        res = "A" + res;
    }
    return res;
}

/**
 * a getter for the Kmer size
 * @return the size of the kmer
 */
ushort Kmer::getSize() {
    return this->size;
}

/**
 * get a subkmer of the original kmer with (0 <= start <= end < k)
 * @param start index of the wanted start of the sub kmer
 * @param end index of the end of the wanted kmer.
 * @return the Kmer[start:end+1]
 */
Kmer Kmer::getSubKmer(int start, int end) {
    uint64_t mask = (1 << ((size - start) * 2)) - 1;
    int decalage = ((size - end - 1) * 2);
    mask = mask >> decalage;
    mask = mask << decalage;
    return Kmer((this->value & mask) >> decalage, (end - start) + 1);
}
