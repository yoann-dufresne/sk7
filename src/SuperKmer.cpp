#include "SuperKmer.hpp"

#include <utility>
#include <iostream>
#include <bitset>

using namespace std;


/**
 * Default constructor for a SuperKmer
 */
SuperKmer::SuperKmer() {
    this->tab = std::vector<TYPE>();
}

/**
 * Constructor for a SuperKmer
 * @param tab the representation of the SuperKmer as a vector of the wanted type
 */
SuperKmer::SuperKmer(std::vector<TYPE> tab) {
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
    ulong section = start / SIZE; // the part of the vector to read in
    while(section >= tab.size()) {
        tab.push_back(0);
    }
    int position = start % SIZE; // the position of the bit to read
    uint64_t mask = 1 << (length - 1);
    for (int i = 0; i < length; i++) {
        uint8_t currentBit = (value & mask) >> (length - 1 - i);
        tab.at(section) &= ~(1 << (SIZE - position - 1)); // reset the old bit
        tab.at(section) |= (currentBit << (SIZE - 1 - position)); // put the new one
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

/**
 * Cut a superKmer in half
 * @param commonPartStart the starting position idx (in nucleotides) of the common part of the cut
 * @param commonPartEnd the ending position idx (in nucleotides) of the common part of the cut
 * @param fixBitSize the size of the prefix and suffix in bits
 * @return a vector containing the two new superKmers
 */
std::vector<SuperKmer> SuperKmer::cut(const int &commonPartStart, const int &commonPartEnd,  const int &fixBitSize) {
    SuperKmer part1 = SuperKmer();
    SuperKmer part2 = SuperKmer();

    int currentPrefixLen = accessBits(0, fixBitSize);
    int currentSuffixLen = accessBits(fixBitSize, 2 * fixBitSize);


    // Shifted to the right
    int part1Suffix = currentSuffixLen - 1;

    part1.setBits(0, fixBitSize, currentPrefixLen);
    part1.setBits(fixBitSize, fixBitSize, part1Suffix);

    uint64_t value1 = 0;
    int cmpt1 = 0;
    while (cmpt1 < part1Suffix || cmpt1 < currentPrefixLen) {
        int readStart = 2 * fixBitSize + cmpt1 * 4;
        if (cmpt1 < part1Suffix) {
            value1 = (value1 << 2) + accessBits(readStart, readStart + 2);
        } else {
            value1 <<= 2;
        }
        readStart = 2 * fixBitSize + cmpt1 * 4 + 2;
        if (cmpt1 < currentPrefixLen) {
            value1 = (value1 << 2) + accessBits(readStart, readStart + 2);
        } else {
            value1 <<= 2;
        }
        cmpt1++;
    }

    part1.setBits(2*fixBitSize, ((part1Suffix < currentPrefixLen)?currentPrefixLen:part1Suffix) * 4, value1);

    // Shifted to the left
    int part2Prefix = currentPrefixLen - 1;

    part2.setBits(0, fixBitSize, part2Prefix);
    part2.setBits(fixBitSize, fixBitSize, currentSuffixLen);

    uint64_t value2 = 0;
    int cmpt2 = 0;
    while (cmpt2 < part2Prefix || cmpt2 < currentSuffixLen) {
        int readStart = 2 * fixBitSize + cmpt2 * 4;
        if (cmpt2 < currentSuffixLen) {
            value2 = (value2 << 2) + accessBits(readStart, readStart + 2);
        } else {
            value2 <<= 2;
        }
        readStart = 2 * fixBitSize + cmpt2 * 4 + 2;
        if (cmpt2 < part2Prefix) {
            value2 = (value2 << 2) + accessBits(readStart, readStart + 2);
        } else {
            value2 <<= 2;
        }
        cmpt2++;
    }

    part2.setBits(2*fixBitSize, ((part2Prefix < currentPrefixLen)?currentSuffixLen:part2Prefix) * 4,value2);

    return std::vector<SuperKmer>({part1, part2});
}

void SuperKmer::print(const int &fixBitSize) {
    cout << "SK print : ";
    for (ulong i = 0; i < tab.size(); i++) {
        cout << "\t" << bitset<8>(tab[i]);
    }
    int prefixLen = getPrefixLen(fixBitSize);
    int suffixLen = getSuffixLen(fixBitSize);
    uint64_t value = getValue(fixBitSize);
    std::string nucleo = Kmer(value, prefixLen + suffixLen).toString();
    if (suffixLen == 0) nucleo = '_' + nucleo;
    int i = 0;
    while(i < (int) nucleo.length()) {
        if(i > 2 * suffixLen) {
            nucleo.at(i) = '_';
        }
        if (i + 1 > 2 * prefixLen) {
            nucleo.at(i + 1) = '_';
        }
        i += 2;
    }

    cout << " --> " << "prefixLen = " << prefixLen << " suffixLen = " << suffixLen << " : "
    << nucleo << endl;

}

/**
 * Access the length of the prefix
 * @param fixBitSize the size (in bits) of the space used to store the prefix and suffix length
 * @return the length of the prefix
 */
int SuperKmer::getPrefixLen(const int &fixBitSize) {
    return accessBits(0, fixBitSize);
}

/**
 * Access the length of the suffix
 * @param fixBitSize the size (in bits) of the space used to store the prefix and suffix length
 * @return the length of the suffix
 */
int SuperKmer::getSuffixLen(const int &fixBitSize) {
    return accessBits(fixBitSize, 2 * fixBitSize);
}

/**
 * Access the value of the SuperKmer
 * @param fixBitSize the size (in bits) of the space used to store the prefix and suffix length
 * @return the value of the SuperKmer
 */
uint64_t SuperKmer::getValue(const int &fixBitSize) {
    int suffixLen = getSuffixLen(fixBitSize);
    int prefixLen = getPrefixLen(fixBitSize);
    return accessBits(2*fixBitSize, 2*fixBitSize + ((suffixLen < prefixLen)?prefixLen:suffixLen) * 4);
}

/**
 * Intersect two superKmers
 * @param toIntersect the SuperKmer to intersect
 * @return a SuperKmer representing the intersection between the two sets of Kmer represented by the two SuperKmers
 */
SuperKmer SuperKmer::operator&(SuperKmer toIntersect) {
    SuperKmer result = SuperKmer();
//    int prefixLen = getPrefixLen();
    return result;
}

