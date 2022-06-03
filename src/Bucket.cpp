#include "src/headers/Bucket.hpp"

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
    int maxLen = (prefixLen < suffixLen)? suffixLen : prefixLen;

    ///START OF THE BINARY SEARCH
    int start = 0;
    int end = this->orderedList.size();
    int fixBitSize = ceil(log2(this->kmerLength - this->minimiserLength + 1));
    int SKheader = 2 * fixBitSize;


    while (start <= end) {
        int middle = (end + start) / 2;
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

        int found = ceil(log2(knownInfo));
        int needed = ceil(log2(withoutMinimiser.getValue()));

        uint64_t toCompare = withoutMinimiser.getValue() & (knownInfo >> (found - needed));

        if (toCompare < maskedSK) {
            end = middle - 1;
            continue;
        }
        if (toCompare > maskedSK) {
            start = middle + 1;
            continue;
        }
        if (toCompare == maskedSK) { //Match
            if (currentPrefixLen >= prefixLen && currentSuffixLen >= suffixLen) {
                return middle;
            } else { //No info
                start += 1;
                continue;
            }
        }
    }

    return -1;
}

