#include "headers/SuperKmer.hpp"

#include <utility>

/**
 * Default constructor for a superKmer
 * @param tab the representation of the SuperKmer as a vector of the wanted type
 */
SuperKmer::SuperKmer(std::vector<uint8_t> tab) {
    this->tab = std::move(tab);
}

/**
 * Return the the bits in [start, end) interval of the superKmer
 * @param start the starting position of the read (included)
 * @param end the ending position of the read (excluded)
 * @return the value of the read bits
 */
uint64_t SuperKmer::accessBits(int start, int end) {
    uint64_t result = 0;
    int reads = end - start; // number if bits to read
    int section = start / SIZE; // the part of the vector to read in
    int position = start % SIZE; // the position of the bit to read
    int mask = 1 << (SIZE - position - 1); // the mask to read the wanted bit
    for (int i = 0; i < reads; i++) {
        result <<= 1;
        result += (this->tab.at(section) & mask) >> (SIZE - position - 1);
        if (++position == SIZE) { // we arrive at the end of the section of the vector
            section++;
            position = 0;
            mask = 1 << (SIZE - position - 1);
        } else { // we continue to read in the same section of the vector
            mask >>= 1;
        }
    }
    return result;
}

