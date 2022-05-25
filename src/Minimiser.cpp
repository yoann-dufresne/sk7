#include "Minimiser.hpp"

/**
 * Constructor for a minimiser
 * @param hashFunction a function that gives an whole value to a kmer
 * @param kmer the kmer to minimise
 * @param size the wanted size of the minimiser
 */
Minimiser ::Minimiser(uint64_t (*hashFunction)(Kmer, ushort), const Kmer& kmer, ushort size) {
    this->hashFunction = hashFunction;
    this->size = size;
    this->value = hashFunction(kmer, size);
}

/**
 * getteur for the minimiser value
 * @return the minimiser's value
 */
uint64_t Minimiser::getValue() const {
    return this->value;
}

std::string Minimiser::toString() const {
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

