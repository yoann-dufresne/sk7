#include <string>
#include <cmath>

#include "headers/Kmer.hpp"


/**
 * Consturctor from value
 * @param value the value of the wanted kmer
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

ushort Kmer::getSize() {
    return this->size;
}
