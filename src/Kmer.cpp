#include <string>
#include <cmath>

#include "Kmer.hpp"

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
