#include "SuperKmer.hpp"

#include <utility>
#include <iostream>

using namespace std;

namespace sk7 {

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
uint64_t SuperKmer::accessBits(int start, int end) const {
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
    while (section >= tab.size()) { // create missing section of the vector
        tab.push_back(0);
    }
    int position = start % SIZE; // the position of the bit to read
    uint64_t mask = (uint64_t) 1 << (length - 1);
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
* Create a mask to represent holes in the SuperKmer (ex : |_C_T -> 00110011)
* @param prefixLen the length of the SuperKmer's prefix
* @param suffixLen the length of the SuperKmer's suffix
* @return the built mask
*/
uint64_t SuperKmer::buildSKMask(const int &prefixLen, const int &suffixLen) const {
    int i = 0;
    uint64_t knownInfo = 0;
    while (i < suffixLen || i < prefixLen) { //Build the mask from superkmer for known information
        if (i < suffixLen) {
            knownInfo = (knownInfo << 2) + 0b11;
        } else {
            knownInfo <<= 2;
        }
        if (i < prefixLen) {
            knownInfo = (knownInfo << 2) + 0b11;
        } else {
            knownInfo <<= 2;
        }
        i++;
    }
    return knownInfo;
}

/**
* Print a SuperKmer
*/
void SuperKmer::print() const {
    if (tab.empty()) {
        cout << "Empty SK" << endl;
        return;
    }
    cout << "SK print : ";
    for (ulong i = 0; i < tab.size(); i++) {
        cout << "\t" << bitset<SIZE>(tab[i]);
    }
    int prefixLen = getPrefixLen();
    int suffixLen = getSuffixLen();
    uint64_t value = getValue();
    std::string nucleo = "";
    for (int i = 0; i < 2 * max(prefixLen, suffixLen); i++) {
        switch (value & 0b11) {
            case 0b00:
                nucleo = "A" + nucleo;
                break;
            case 0b01:
                nucleo = "C" + nucleo;
                break;
            case 0b10:
                nucleo = "T" + nucleo;
                break;
            case 0b11:
                nucleo = "G" + nucleo;
                break;
        }
        value >>= 2;
    }
    for (int i = 0; i < 2 * max(prefixLen, suffixLen); i++) { // replace excess nucleotides by "_"
        if (i + 1 > 2 * suffixLen) {
            nucleo.at(i) = '_';
        }
        i++;
        if (i + 1 > 2 * prefixLen) {
            nucleo.at(i) = '_';
        }
    }
    cout << " --> " << "prefixLen = " << prefixLen << ", suffixLen = " << suffixLen << " : "
         << nucleo << endl;

}

/**
* Access the length of the prefix
* @return the length of the prefix
*/
int SuperKmer::getPrefixLen() const {
    return accessBits(0, sk7::fixBitSize);
}

/**
* Access the length of the suffix
* @return the length of the suffix
*/
int SuperKmer::getSuffixLen() const {
    return accessBits(sk7::fixBitSize, 2 * sk7::fixBitSize);
}

/**
* Access the value of the SuperKmer
* @return the value of the SuperKmer
*/
uint64_t SuperKmer::getValue() const {
    int suffixLen = getSuffixLen();
    int prefixLen = getPrefixLen();
    return accessBits(2 * sk7::fixBitSize, 2 * sk7::fixBitSize + max(prefixLen, suffixLen) * 4);
}

/**
* Check if two SuperKmers are equal
* @param toCompare the SuperKmer to compare
* @return true if the two SuperKmer are equal, else false
*/
bool SuperKmer::operator==(const SuperKmer &toCompare) const {
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
    // getting needed info
    SuperKmer result = SuperKmer();
    int prefixLen = getPrefixLen();
    int suffixLen = getSuffixLen();
    int maxLen = max(prefixLen, suffixLen);
    uint64_t value = getValue();

    int prefixLenIntersect = toIntersect.getPrefixLen();
    int suffixLenIntersect = toIntersect.getSuffixLen();
    int maxLenIntersect = max(prefixLenIntersect, suffixLenIntersect);
    uint64_t valueIntersect = toIntersect.getValue();

    if (maxLen < maxLenIntersect) { // aligning the values
        valueIntersect >>= 4 * (maxLenIntersect - maxLen);
    } else {
        value >>= 4 * (maxLen - maxLenIntersect);
    }

    // Separation of the prefix and suffix
    uint64_t prefixMask = 0b00;
    uint64_t suffixMask = 0b00;
    for (int i = 0; i < 2 * min(maxLen, maxLenIntersect); i++) {
        if (i % 2 == 0) {
            prefixMask <<= 2;
            suffixMask = (suffixMask << 2) + 0b11;
        } else {
            prefixMask = (prefixMask << 2) + 0b11;
            suffixMask <<= 2;
        }
    }

    uint64_t prefixValue = value & prefixMask;
    uint64_t suffixValue = value & suffixMask;
    uint64_t prefixValueIntersect = valueIntersect & prefixMask;
    uint64_t suffixValueIntersect = valueIntersect & suffixMask;

    uint64_t prefixComparison = prefixValue ^ prefixValueIntersect;
    uint64_t suffixComparison = suffixValue ^ suffixValueIntersect;

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

    // searching for the first difference in suffix and prefix
    int firstDiffIndexPrefix = leftmost1IndexPrefix - startOfMeaningfulIndex;
    int firstDiffNucleotidePrefix = min({firstDiffIndexPrefix / 4, prefixLen, prefixLenIntersect});

    int firstDiffIndexSuffix = leftmost1IndexSuffix - startOfMeaningfulIndex;
    int firstDiffNucleotideSuffix = min({firstDiffIndexSuffix / 4, suffixLen, suffixLenIntersect});

    // Empty intersection
    if (firstDiffNucleotideSuffix + firstDiffNucleotidePrefix + sk7::m < sk7::k) {
        return SuperKmer();
    }

    // building intersection
    result.setBits(0, sk7::fixBitSize, firstDiffNucleotidePrefix);
    result.setBits(sk7::fixBitSize, sk7::fixBitSize, firstDiffNucleotideSuffix);

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
    // Only works for small SK due to uint64_t
    uint64_t final_value = (((value & prefixMask) + (value & suffixMask)) >>
                                                                          (4 * (min(maxLen, maxLenIntersect) -
                                                                                max(firstDiffNucleotidePrefix,
                                                                                    firstDiffNucleotideSuffix))));
    result.setBits(2 * sk7::fixBitSize, 4 * max(firstDiffNucleotidePrefix, firstDiffNucleotideSuffix),
                   final_value);

    return result;
}

/**
* Unite two compatible SuperKmers
* @param toUnite the SuperKmer to unite
* @return the union of the two SuperKmers
*/
SuperKmer SuperKmer::operator|(SuperKmer toUnite) {
    SuperKmer result = SuperKmer();
    int newPrefixLen = max(getPrefixLen(), toUnite.getPrefixLen());
    result.setBits(0, sk7::fixBitSize, newPrefixLen);
    int newSuffixLen = max(getSuffixLen(), toUnite.getSuffixLen());
    result.setBits(sk7::fixBitSize, sk7::fixBitSize, newSuffixLen);

    uint64_t value;
    int diff = max(getPrefixLen(), getSuffixLen()) -
               max(toUnite.getSuffixLen(), toUnite.getPrefixLen()); // difference in size (in bits)
    if (diff > 0) {
        value = getValue() | (toUnite.getValue() << 4 * diff);
    } else {
        value = (getValue() << 4 * (-diff)) | (toUnite.getValue());
    }

    result.setBits(2 * sk7::fixBitSize, 4 * max(newPrefixLen, newSuffixLen), value);
    return result;
}


/**
* Translate a SuperKmer to a list of SuperKmers that represent each a single Kmer
* @return the list of interleaved Kmers represented by the SuperKmer
*/
std::vector<SuperKmer> SuperKmer::split() const {

    int prefixLen = getPrefixLen();
    int suffixLen = getSuffixLen();

    if (prefixLen + sk7::m > sk7::k || suffixLen + sk7::m > sk7::k || prefixLen + suffixLen + sk7::m < sk7::k) {
        throw runtime_error("Incoherent SuperKmer with set k and m");
    }

    std::vector<SuperKmer> result = vector<SuperKmer>(sk7::k - sk7::m + 1);

    int maxPref = sk7::k - sk7::m;
    int nbKmer = sk7::k - sk7::m + 1;

    for (int i = 0; i < nbKmer; i++) {

        int currentPrefixLen = maxPref - i;
        int currentSuffixLen = i;

        if (currentPrefixLen > prefixLen || currentSuffixLen > suffixLen) {
            result.at(i) = SuperKmer();
            continue;
        }
        SuperKmer toInsert = SuperKmer();

        Kmer currentKmer = readKmer(currentPrefixLen);
        uint64_t sum = currentKmer.getValue();

        toInsert.setBits(0, sk7::fixBitSize, currentPrefixLen);
        toInsert.setBits(sk7::fixBitSize, sk7::fixBitSize, currentSuffixLen);
        toInsert.setBits(sk7::fixBitSize * 2, 4 * max(currentSuffixLen, currentPrefixLen), sum);
        result.at(currentSuffixLen) = toInsert;
    }
    return result;
}

/**
* Compare every Kmers included in two SuperKmers
* @deprecated not up to date
* @param superKmer1 the first SuperKmer
* @param superKmer2 the second SuperKmer
* @return a vector with the result of every comparison in order
*/
vector <SuperKmer::logic> SuperKmer::compareSK(SuperKmer &superKmer1, SuperKmer &superKmer2) {
    std::vector<SuperKmer::logic> result = std::vector<SuperKmer::logic>(sk7::k - sk7::m + 1);

    std::vector<SuperKmer> unpacked1 = superKmer1.split();
    std::vector<SuperKmer> unpacked2 = superKmer2.split();

    for (int i = 0; i < sk7::k - sk7::m + 1; ++i) {
        int prefix1Len = unpacked1.at(i).getPrefixLen();
        int prefix2Len = unpacked2.at(i).getPrefixLen();

        if ((!prefix1Len && !unpacked1.at(i).getSuffixLen()) ||
            (!prefix2Len && !unpacked2.at(i).getSuffixLen())) { // One is empty
            result.at(i) = SuperKmer::INCOMPARABLE;
            continue;
        }

        if (prefix1Len < prefix2Len) { // Not the same minimiser position
            result.at(i) = SuperKmer::INCOMPARABLE;
            continue;
        }
        if (prefix1Len > prefix2Len) {
            result.at(i) = SuperKmer::INCOMPARABLE;
            continue;

        } else { // Same minimiser position

            uint64_t value1 = unpacked1.at(i).getValue();
            uint64_t value2 = unpacked2.at(i).getValue();

            if (value1 == value2) {
                result.at(i) = SuperKmer::EQUAL;
            } else if (value1 > value2) {
                result.at(i) = SuperKmer::SUPERIOR;
            } else {
                result.at(i) = SuperKmer::INFERIOR;
            }
            continue;
        }
    }

    return result;
}

/**
* Find the value of the non-interleaved Kmer of a SuperKmer that contains an unique Kmer
* @return the value of the non-interleaved Kmer found with no minimiser
*/
uint64_t SuperKmer::nonInterleavedKmerValue() const {
    int prefixLen = getPrefixLen(); // superkmer's prefix length
    int suffixLen = getSuffixLen(); // superkmer's suffix length
    uint64_t kmerValue = 0;
    for (int i = prefixLen - 1; i >= 0; i--) {
        int readStart = 2 * sk7::fixBitSize + i * 4 + 2;
        kmerValue = (kmerValue << 2) + accessBits(readStart, readStart + 2);
    }
//    kmerValue = (kmerValue << sk7::m * 2); // make space for the minimiser
    for (int i = 0; i < suffixLen; i++) {
        int readStart = 2 * sk7::fixBitSize + i * 4;
        kmerValue = (kmerValue << 2) + accessBits(readStart, readStart + 2);
    }
    return kmerValue;
}

/**
* Read a specific Kmer from the SuperKmer
* @param kmerPrefixLen the length of the prefix of the wanted Kmer
* @return the value of the Kmer (stripped from minimiser and interleaved)
*/
Kmer SuperKmer::readKmer(int kmerPrefixLen) const {
    int prefixLen = getPrefixLen();
    int suffixLen = getSuffixLen();

    // Not a corresponding Kmer
    if (prefixLen < kmerPrefixLen || kmerPrefixLen > sk7::k - sk7::m ||
        kmerPrefixLen + suffixLen < sk7::k - sk7::m) {
        return Kmer(0, 0);
    }

    int maxLen = max(prefixLen, suffixLen);
    uint64_t value = getValue();

    // Separation of the prefix and suffix with a mask
    uint64_t prefixMask = 0b00;
    uint64_t suffixMask = 0b00;
    for (int i = 0; i < 2 * maxLen; i++) {
        if (i % 2 == 0) {
            prefixMask <<= 2;
            suffixMask = (suffixMask << 2) + 0b11;
        } else {
            prefixMask = (prefixMask << 2) + 0b11;
            suffixMask <<= 2;
        }
    }

    suffixMask = suffixMask << 4 * (maxLen - (sk7::k - sk7::m - kmerPrefixLen));
    prefixMask = prefixMask << 4 * (maxLen - kmerPrefixLen);

    uint64_t forPref = value & prefixMask;
    uint64_t forSuf = value & suffixMask;

    // Build the result
    uint64_t sum = (forSuf + forPref) >> 4 * (maxLen - max(kmerPrefixLen, sk7::k - sk7::m - kmerPrefixLen));
    return Kmer(sum, sk7::k - sk7::m);
}

/**
* Extract a SuperKmer that represents a single Kmer from this SuperKmer
* @param wantedPrefixLen
* @return
*/
SuperKmer SuperKmer::extract(int wantedPrefixLen) const {
    int prefixLen = getPrefixLen();
    int suffixLen = getSuffixLen();

    if (prefixLen < wantedPrefixLen || wantedPrefixLen > sk7::k - sk7::m ||
        wantedPrefixLen + suffixLen < sk7::k - sk7::m) {
        return SuperKmer();
    }

    SuperKmer result = SuperKmer();

    int maxLen = max(prefixLen, suffixLen);
    uint64_t value = getValue();

    // Separation of the prefix and suffix with a mask
    uint64_t prefixMask = 0b00;
    uint64_t suffixMask = 0b00;
    for (int i = 0; i < 2 * maxLen; i++) {
        if (i % 2 == 0) {
            prefixMask <<= 2;
            suffixMask = (suffixMask << 2) + 0b11;
        } else {
            prefixMask = (prefixMask << 2) + 0b11;
            suffixMask <<= 2;
        }
    }

    suffixMask = suffixMask << 4 * (maxLen - (sk7::k - sk7::m - wantedPrefixLen));
    prefixMask = prefixMask << 4 * (maxLen - wantedPrefixLen);

    uint64_t forPref = value & prefixMask;
    uint64_t forSuf = value & suffixMask;

    uint64_t sum = (forSuf + forPref) >> 4 * (maxLen - max(wantedPrefixLen, sk7::k - sk7::m - wantedPrefixLen));

    result.setBits(0, sk7::fixBitSize, wantedPrefixLen);
    result.setBits(sk7::fixBitSize, sk7::fixBitSize, sk7::k - sk7::m - wantedPrefixLen);
    result.setBits(2 * sk7::fixBitSize, 4 * max(wantedPrefixLen, sk7::k - sk7::m - wantedPrefixLen), sum);

    return result;
}

}