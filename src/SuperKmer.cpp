#include "SuperKmer.hpp"

#include <utility>
#include <iostream>
#include <bitset>

using namespace std;

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



/**
 * Set the bits in the range [start, start + length) with those of value in the same order
 * @param start the position of the starting bit
 * @param length the length of the modification
 * @param value the value to change the bits with
 */
void SuperKmer::setBits(const int &start, const int &length, const uint64_t &value) {
    int section = start / SIZE; // the part of the vector to read in
    int position = start % SIZE; // the position of the bit to read
    uint64_t mask = 1 << (length - 1);
    for (int i = 0; i < length; i++) {
        uint8_t currentBit = (value & mask) >> (length - 1 - i);
        tab.at(section) &= ~(1 << (SIZE - position - 1));
        tab.at(section) |= (currentBit << (SIZE - 1 - position));
        mask >>= 1;
        if (++position == SIZE) { // we arrive at the end of the section of the vector
            if ((ulong) ++section == tab.size()) { // Make a new section of the tab
                tab.push_back(0);
            }
            position = 0;
        } else { // we continue to read in the same section of the vector
            continue;
        }
    }
}

void SuperKmer::print() {
    cout << "SK print : ";
    for (ulong i = 0; i < tab.size(); i++) {
        cout << "\t" << bitset<8>(tab[i]);
    }
    cout << endl;
}

