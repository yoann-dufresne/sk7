#include "headers/Minimiser.hpp"

/**
 * Constructor for a minimiser
 * @param hashFunction a function that gives an whole value to a mmer
 * @param size the wanted size of the minimiser
 */
Minimiser ::Minimiser(hashPos (*hashFunction)(Kmer, ushort), ushort size) {
    this->hashFunction = hashFunction;
    this->size = size;
    this->value = 0;
    this->pos = 0;
}

/**
 * Initialise value with the hash of the given kmer
 * @param kmer the initializer
 */
void Minimiser::init(Kmer kmer) {
    hashPos retHash = hashFunction(kmer, this->size);
    this->value = retHash.hashValue;
    this->pos = retHash.pos;
}

/**
 * Build a minimiser depending on the new end
 * @param hashFunction a function that gives an whole value to a mmer
 * @param kmer the current kmer to minimise
 * @param size the wanted size of the minimiser
 * @return the current minimiser
 */
void Minimiser::fromNewEnd(Kmer kmer) {
    if (pos == 0) {
        this->init(kmer);
        return;
    }
    Kmer end = kmer.getSubKmer(kmer.getSize() - size, kmer.getSize() - 1);
    hashPos retHash = hashFunction(end, size);
    uint64_t end_val = retHash.hashValue;
    if (end_val < this->value) { //The new ending is the minimiser
        this->value = end_val;
        this->pos = kmer.getSize() - size;
    } else {
        this->pos--; //We keep the same minimiser but its position shift
    }
}

/**
 * Getter for the minimiser value
 * @return the minimiser's value
 */
uint64_t Minimiser::getValue() const {
    return this->value;
}

/**
 * Getter for the pos attribute
 * @return the minimiser's pos
 */
short Minimiser::getPos() const {
    return this->pos;
}

/**
 * Give a textual representation of a minimiser
 * @return a string on {'A', 'C', 'G', 'T'} that represents the minimiser
 */
std::string Minimiser::toString() const {
    uint64_t tmp = value;
    std::string res;
    res.reserve(this->size);
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
    while (res.length() < this->size) {
        res = "A" + res;
    }
    return res;
}

