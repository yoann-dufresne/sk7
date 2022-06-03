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
 * @return true if kmer is in the Bucket else false
 */
bool Bucket::isIn(Kmer kmer) {
    ///PREPARATION OF THE SEARCH
    cout << endl << endl;
//    if (kmer.getLength() < this->kmerLength) {
//        throw (std::runtime_error(std::string("Illegal length for the given Kmer")));
//    }
    Minimiser kmerMinimiser = Minimiser(alpha, this->minimiserLength, kmer);
    Kmer withoutMinimiser = kmer.removePart(kmerMinimiser.getPos(), this->minimiserLength);
    cout << "striped : " << withoutMinimiser.toString() << endl;

    uint64_t kmerMask = interleavedOrder(withoutMinimiser, kmerMinimiser.getPos());

    cout << "ordered : " << withoutMinimiser.toString() << endl;
    cout << "kmerMask : " << kmerMask << endl;
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
//        int middle = 2; //For test
        SuperKmer current = this->orderedList.at(middle); //Current SuperKmer of the list

        int currentPrefixLen = current.accessBits(0, fixBitSize); //Current superkmer's prefix length
        int currentSuffixLen = current.accessBits(fixBitSize, SKheader); //Current superkmer's suffix length
        int currentMaxLen = (currentPrefixLen < currentSuffixLen) ? currentSuffixLen : currentPrefixLen;

        cout << " Current prefix length : " << currentPrefixLen << " current suffix length : " << currentSuffixLen
             << endl;


        int totalNucleotides = currentPrefixLen + currentSuffixLen; //Taille du superkmer
        uint64_t currentValue = current.accessBits(SKheader, SKheader + currentMaxLen * 4); //Valeur du superkmer

        cout << "current value : " << currentValue << endl;

        uint64_t maskedSK = (currentValue >> ((currentMaxLen - maxLen) * 4)) & kmerMask;
        cout << Kmer(currentValue, totalNucleotides).toString() << " kmer masked to " << maskedSK << endl;

        uint64_t knownInfo = 0;
        int i = 0;
        int needed = ceil(log2(withoutMinimiser.getValue()));
        while (true) {
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

        cout << "knownInfo : " << knownInfo << endl;
        int found = ceil(log2(knownInfo));

        cout << withoutMinimiser.getValue() << endl;
        cout << "found : " << found << " needed : " << needed << endl;

        uint64_t toCompare = withoutMinimiser.getValue() & (knownInfo >> (found - needed));
        cout << "To compare : " << toCompare << " with : " << maskedSK << endl;

        if (toCompare < maskedSK) {
            end = middle - 1;
            continue;
        }
        if (toCompare > maskedSK) {
            start = middle + 1;
            continue;
        }
        if (toCompare == maskedSK) {
            if (currentPrefixLen >= prefixLen && currentSuffixLen >= suffixLen) {
                return true;
            } else {
                start += 1;
                continue;
            }
        }

    }

    return false;
}

