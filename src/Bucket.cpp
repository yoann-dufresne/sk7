#include "headers/Bucket.hpp"

#include <iostream>
#include <cmath>

using namespace std;

/**
 * Initialise a Bucket
 * @param minimiserLength the length of the minimiser linked to the bucket
 * @param minimiserValue the value of the minimiser
 * @param orderedList the list of superKmer in the bucket
 * @param kmerLength the length of the kmer used to create the superkmer
 */
Bucket::Bucket(int minimiserLength, uint64_t minimiserValue, int kmerLength) {
    this->minimiserLength = minimiserLength;
    this->minimiser = minimiserValue;
    this->orderedList = std::vector<SuperKmer>();
    this->kmerLength = kmerLength;
}

/**
 * Add a element to the Bucket (for now doesn't consider order)
 * @param superKmer the superkmer to add to the bucket
 */
void Bucket::addToList(SuperKmer superKmer) {
    this->orderedList.push_back(superKmer);
}


/**
 * Search a element in a bucket
 * @param kmer the Kmer to search
 * @return the position of the kmer in the list if found else -1
 */
int Bucket::isIn(Kmer kmer) {
    ///PREPARATION OF THE SEARCH
    Minimiser kmerMinimiser = Minimiser(alpha, this->minimiserLength, kmer);
    Kmer withoutMinimiser = kmer.removePart(kmerMinimiser.getPos(), this->minimiserLength);

    uint64_t kmerMask = interleavedOrder(withoutMinimiser, kmerMinimiser.getPos());

    int prefixLen = kmerMinimiser.getPos();
    int suffixLen = kmer.getLength() - kmerMinimiser.getPos() - this->minimiserLength;
    int maxLen;

    int needed; // number of bits to compare
    if (prefixLen > suffixLen) {
        maxLen = prefixLen;
        needed = maxLen * 4;
    } else {
        maxLen = suffixLen;
        needed = maxLen * 4 - 2;
    }
    if (prefixLen == suffixLen) {
        needed = 4 * prefixLen;
    }

    int fixBitSize = ceil(log2(this->kmerLength - this->minimiserLength + 1));
    int SKheader = 2 * fixBitSize;

    ///START OF THE BINARY SEARCH
    int start = 0;
    int end = this->orderedList.size() - 1;
    bool infoFound = true;
    int lastPositionWithInformation;
    int lastStartWithInformation = 0;
    while (start <= end) {
        int middle = (end + start) / 2;
        if (infoFound) {
            lastPositionWithInformation = middle;
            infoFound = false;
        }
        SuperKmer current = this->orderedList.at(middle); //Current SuperKmer of the list

        int currentPrefixLen = current.accessBits(0, fixBitSize); //Current superkmer's prefix length
        int currentSuffixLen = current.accessBits(fixBitSize, SKheader); //Current superkmer's suffix length
        int currentMaxLen = (currentPrefixLen < currentSuffixLen) ? currentSuffixLen : currentPrefixLen;

        uint64_t currentValue = current.accessBits(SKheader, SKheader + currentMaxLen * 4); //Valeur du superkmer

        uint64_t maskedSK = (currentValue >> ((currentMaxLen - maxLen) * 4)) & kmerMask;

        uint64_t knownInfo = 0;
        int i = 0;

        while (true) { //Build the mask from superkmer for known information
            if(i < currentSuffixLen) {
                knownInfo = (knownInfo << 2) + 0b11;
            } else {
                knownInfo <<=2;
            }
            if (i < currentPrefixLen) {
                knownInfo = (knownInfo << 2) + 0b11;
            } else {
                knownInfo <<=2;
            }
            i++;
            if (i >= currentSuffixLen && i >= currentPrefixLen) break;
        }

        int found = 4 * i; //ceil(log2(knownInfo));

        //We align the values
        uint64_t toCompare;
//        cout << "middle : " << middle << " found : " << found << " needed : " << needed << " info : " << knownInfo << endl;
//        cout << "before : " << withoutMinimiser.getValue() << " " <<((knownInfo >> (found - needed))) << endl;
//        cout << "string : " << withoutMinimiser.toString() << endl;
        if (found >= needed) {
            if (prefixLen < suffixLen) toCompare = withoutMinimiser.getValue() & (knownInfo >> (found - needed - 2));
            else toCompare = withoutMinimiser.getValue() & (knownInfo >> (found - needed));
        } else {
            toCompare = withoutMinimiser.getValue() & (knownInfo << (needed - found));
        }

//        cout << "We're comparing : " << toCompare << " and : " << maskedSK << endl;
        if (toCompare < maskedSK) {
            end = lastPositionWithInformation - 1;
            infoFound = true;
            continue;
        }

        if (toCompare > maskedSK) {
            start = lastPositionWithInformation + 1;
            lastStartWithInformation = start;
            infoFound = true;
            continue;
        }

        if (toCompare == maskedSK) { //Match or not enough information
            if (currentPrefixLen >= prefixLen && currentSuffixLen >= suffixLen) { // Match
                return middle;
            } else { //No info
                if (start < end) {
                    start += 1;
                } else { // No info till the end of the linear search
                    start = lastStartWithInformation;
                    end = lastPositionWithInformation - 1;
                    infoFound = true;
                }
                continue;
            }
        }
    }

    return -1;
}
