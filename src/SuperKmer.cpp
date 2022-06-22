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
    tab.push_back(0);
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
    while(section >= tab.size()) { // create missing section of the vector
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
            if ((ulong) ++section == tab.size() && i != length - 1) { // Make a new section of the tab
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

void SuperKmer::print() {
    if (tab.empty()) {
        cout << "Empty SK" << endl;
        return;
    }
    cout << "SK print : ";
    for (ulong i = 0; i < tab.size(); i++) {
        cout << "\t" << bitset<SIZE>(tab[i]);
    }
    int prefixLen = getPrefixLen(sk7::fixBitSize);
    int suffixLen = getSuffixLen(sk7::fixBitSize);
    uint64_t value = getValue(sk7::fixBitSize);
    std::string nucleo = Kmer(value, prefixLen + suffixLen).toString();
    if (suffixLen == 0) nucleo = '_' + nucleo;
    int i = 0;
    while(i < (int) nucleo.length()) {
        if(i >= 2 * suffixLen) {
            nucleo.at(i) = '_';
        }
        try {
            if (i + 1 >= 2 * prefixLen) {
                nucleo.at(i + 1) = '_';
            }
        } catch (...){
            nucleo += '_';
        }
        i += 2;
    }

    cout << " --> " << "prefixLen = " << prefixLen << ", suffixLen = " << suffixLen << " : "
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
 * Check if two SuperKmers are equal
 * @param toCompare the SuperKmer to compare
 * @return true if the two SuperKmer are equal, else false
 */
bool SuperKmer::operator==(const SuperKmer toCompare) {
    if (tab.size() != toCompare.tab.size()) {
        return false;
    }
    for (uint64_t i = 0; i < tab.size(); i++) {
        if (tab.at(i) != toCompare.tab.at(i)) {
            return false;
        }
    }
    return true;
}



/**
 * Intersect two superKmers
 * @param toIntersect the SuperKmer to intersect
 * @return a SuperKmer representing the intersection between the two sets of Kmer represented by the two SuperKmers
 */
SuperKmer SuperKmer::operator&(SuperKmer toIntersect) {
    SuperKmer result = SuperKmer();
    int prefixLen = getPrefixLen(sk7::fixBitSize);
    int suffixLen = getSuffixLen(sk7::fixBitSize);
    int maxLen = max(prefixLen, suffixLen);
    uint64_t value = getValue(sk7::fixBitSize);
//    cout << "value : " << value << endl;

    int prefixLenIntersect = toIntersect.getPrefixLen(sk7::fixBitSize);
    int suffixLenIntersect = toIntersect.getSuffixLen(sk7::fixBitSize);
    int maxLenIntersect = max(prefixLenIntersect, suffixLenIntersect);
    uint64_t valueIntersect = toIntersect.getValue(sk7::fixBitSize);
//    cout << "valueIntersect = " << valueIntersect << endl;

    if (maxLen < maxLenIntersect) { // aligning the values
        valueIntersect >>= 4 * (maxLenIntersect - maxLen);
    } else {
        value >>= 4 * (maxLen - maxLenIntersect);
    }

    // Separation of the prefix and suffix
    uint64_t prefixMask = 0b00;
    uint64_t suffixMask = 0b00;
    for (int i = 0; i < 2*min(maxLen, maxLenIntersect); i++) {
        if (i % 2 == 0) {
            prefixMask <<= 2;
            suffixMask = (suffixMask << 2) + 0b11;
        } else {
            prefixMask = (prefixMask << 2) + 0b11;
            suffixMask <<= 2;
        }
    }

//    cout << "prefix mask : " << prefixMask << " suffix mask : " << suffixMask << endl;
    uint64_t prefixValue = value & prefixMask;
    uint64_t suffixValue = value & suffixMask;
    uint64_t prefixValueIntersect = valueIntersect & prefixMask;
    uint64_t suffixValueIntersect = valueIntersect & suffixMask;

//    cout << "prefix val : " << prefixValue << " suffix val : " << suffixValue << endl;
//    cout << "prefix val inter : " << prefixValueIntersect << " suffix val inter : " << suffixValueIntersect <<endl;

    uint64_t prefixComparison = prefixValue ^ prefixValueIntersect;
    uint64_t suffixComparison = suffixValue ^ suffixValueIntersect;
//    cout << "pref comp : " << prefixComparison << endl;
//    cout << "suff comp : " << suffixComparison << endl;
    // Find the leftmost different nucleotide
    int shiftsPrefix = 0;
    while (prefixComparison != 0) { // Finding leftmost 1 in the xor
        prefixComparison >>= 1;
        shiftsPrefix++;
    }

    int shiftsSuffix = 0;
    while (suffixComparison != 0) { // Finding leftmost 1 in the xor
        suffixComparison >>= 1;
        shiftsSuffix++;
    }

    int leftmost1IndexPrefix = 64 - shiftsPrefix;
    int leftmost1IndexSuffix = 64 - shiftsSuffix;

//    cout << "left 1 pref = " << leftmost1IndexPrefix << " left 1 suff : " << leftmost1IndexSuffix << endl;

    int startOfMeaningfulIndex = 64 - min(maxLen, maxLenIntersect) * 4;

    int firstDiffIndexPrefix = leftmost1IndexPrefix - startOfMeaningfulIndex;
    int firstDiffNucleotidePrefix = min({firstDiffIndexPrefix / 4,  prefixLen,  prefixLenIntersect});

    int firstDiffIndexSuffix = leftmost1IndexSuffix - startOfMeaningfulIndex;
    int firstDiffNucleotideSuffix = min({firstDiffIndexSuffix / 4,  suffixLen,  suffixLenIntersect});

    result.setBits(0, sk7::fixBitSize, firstDiffNucleotidePrefix);
    result.setBits(sk7::fixBitSize,  sk7::fixBitSize, firstDiffNucleotideSuffix);

//    cout << "1st diff pref = " << firstDiffNucleotidePrefix << " suff : " << firstDiffNucleotideSuffix << endl;

    prefixMask = 0;
    suffixMask = 0;
    for (int i = 0; i < min(maxLen, maxLenIntersect); i++) {
        if (i >= firstDiffNucleotideSuffix) {
            suffixMask <<= 4;
        } else {
            suffixMask = (suffixMask << 4) + 0b1100;
        }
        if (i >= firstDiffNucleotidePrefix) {
            prefixMask <<= 4;
        } else {
            prefixMask = (prefixMask << 4) + 0b1100;
        }
    }
    prefixMask >>= 2;
//    cout << "new mask pref : " << prefixMask << " suffix : " << suffixMask << endl;
    uint64_t final_value = (((value & prefixMask) + (value & suffixMask)) >> (4*(maxLen - max(firstDiffNucleotidePrefix,  firstDiffNucleotideSuffix))));
//    cout << "f value : " << final_value << " size : " << 4*max(firstDiffNucleotidePrefix, firstDiffNucleotideSuffix) << endl;
    result.setBits(2*sk7::fixBitSize, 4*max(firstDiffNucleotidePrefix, firstDiffNucleotideSuffix),
                   final_value);

//    result.print();
    return result;
}

