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
 * Add a element to the Bucket WITHOUT considering the order
 * @param superKmer the superkmer to add to the bucket
 */
void Bucket::addToList(SuperKmer superKmer) {
    this->orderedList.push_back(superKmer);
}

/**
 * Binary search a kmer in the bucket to find it's position or a viable one
 * @param kmer the kmer to search
 * @param position a parameter used to store the position found by the search
 * @return true if the kmer was already in the bucket false otherwise
 */
bool Bucket::find(Kmer kmer, int &position) {
//    cout << endl << endl;
//    cout << "searching : " << kmer.toString()  << endl;
    ///PREPARATION OF THE SEARCH
    Minimiser kmerMinimiser = Minimiser(alpha, this->minimiserLength, kmer);
    Kmer withoutMinimiser = kmer.removePart(kmerMinimiser.getPos(), this->minimiserLength);

    uint64_t kmerMask = interleavedOrder(withoutMinimiser, kmerMinimiser.getPos());
//    cout << "striped & interleaved : " << withoutMinimiser.toString() << " of value : " << withoutMinimiser.getValue() << endl;
//    cout << "mask = " << kmerMask << endl;

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
    bool lastWasSuperior = false;
    while (start <= end) {
        int middle = (end + start) / 2;
//        cout << "start : " << start << " end : " << end << " middle : " << middle << endl;
        if (infoFound) {
            lastPositionWithInformation = middle;
            infoFound = false;
        }
        SuperKmer current = this->orderedList.at(middle); //Current SuperKmer of the list

        int currentPrefixLen = current.accessBits(0, fixBitSize); //Current superkmer's prefix length
        int currentSuffixLen = current.accessBits(fixBitSize, SKheader); //Current superkmer's suffix length
        int currentMaxLen = (currentPrefixLen < currentSuffixLen) ? currentSuffixLen : currentPrefixLen;

        uint64_t currentValue = current.accessBits(SKheader, SKheader + currentMaxLen * 4); // superkmer's value
//        //if (currentPrefixLen < currentSuffixLen) currentValue >>=2;
//        cout << "current sk : " << Kmer(currentValue, currentSuffixLen + currentPrefixLen).toString() << " of value : " << currentValue <<endl;
        uint64_t maskedSK;
//        cout << "maxLen : " << maxLen << " current maxLen : " << currentMaxLen << endl;

        if (currentMaxLen < maxLen) {
            maskedSK = currentValue & (kmerMask >> ((maxLen - currentMaxLen) * 4));
//            cout << "effective mask = " << (kmerMask >> ((maxLen - currentMaxLen) * 4)) << endl;
        }
        else {
//            cout << "effective comparator = " << (currentValue >> ((currentMaxLen - maxLen) * 4)) << endl;
            maskedSK = (currentValue >> ((currentMaxLen - maxLen) * 4)) & kmerMask;
        }

        uint64_t knownInfo = 0;
        int i = 0;

//        cout << "current SuffixLen : " << currentSuffixLen << " current PrefixLen : " << currentPrefixLen << endl;
        while (i < currentSuffixLen || i < currentPrefixLen) { //Build the mask from superkmer for known information
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
        }

        int found;
        if (currentPrefixLen < currentSuffixLen) {
            found = 4 * i - 2;
            knownInfo >>= 2;
//            maskedSK >>= 2; //could be used to remove excess 0b00
        } else {
            found = 4 * i;
        }
//        cout << "known info : " << knownInfo << endl;


        //We align the values
//        cout << "found : " << found << " needed : " << needed << endl;
//        cout << "prefixLen : " << prefixLen << " suffixLen : " << suffixLen << endl;
        uint64_t toCompare;
        if (found == needed) {
            if (prefixLen < suffixLen) {
                toCompare = (withoutMinimiser.getValue() >> 2) & knownInfo;
                maskedSK >>=2;
            }
            else toCompare = withoutMinimiser.getValue() & knownInfo;
        }
        else if (found > needed) {
            if (prefixLen < suffixLen) toCompare = withoutMinimiser.getValue() & (knownInfo >> (found - needed - 2));
            else toCompare = withoutMinimiser.getValue() & (knownInfo >> ((found - needed - 2)));
//            /*maskedSK >>= (found - needed);*/
        } else { // (found < needed)
            if (prefixLen < suffixLen) {
//                cout << "init value : " << withoutMinimiser.getValue() << " shifted : " << (withoutMinimiser.getValue() >> (needed - found)) << endl;
                toCompare = (withoutMinimiser.getValue() >> (needed - found)) & knownInfo;
            }
            else {
//                cout << "init value : " << withoutMinimiser.getValue() << " shifted : " << (withoutMinimiser.getValue() >> (needed - found)) << endl;
                toCompare = (withoutMinimiser.getValue() >> (needed - found)) & knownInfo;
            }
            if (currentPrefixLen < currentSuffixLen) maskedSK >>=2;

        }

//        cout << "toCompare : " << toCompare << " maskedSK : " << maskedSK << endl;
        if (toCompare < maskedSK) {
//            cout << "before " << endl << endl;
            start = lastStartWithInformation;
            end = lastPositionWithInformation - 1;
            infoFound = true;
            lastWasSuperior = true;
            continue;
        }

        if (toCompare > maskedSK) {
//            cout << "after " << endl << endl;
            start = lastPositionWithInformation + 1;
            lastStartWithInformation = start;
            infoFound = true;
            lastWasSuperior = false;
            continue;
        }

        if (toCompare == maskedSK) { //Match or not enough information
            if (currentPrefixLen >= prefixLen && currentSuffixLen >= suffixLen) { // Match
                position = middle;
                return true;
            } else { //No info
//                cout << "no info" << endl << endl;
                if (start < end) {
                    start += 1;
                    continue;
                }
                if (start == end) { // End of search via linear search
                    position = start;
                    return false;
                }
                else { // No info till the end of the linear search
                    start = lastStartWithInformation;
                    end = lastPositionWithInformation - 1;
                    infoFound = true;
                }
                continue;
            }
        }
    }

//    cout << "ending start : " << start << " ending end : " << end << endl;
//    cout << "last was sup : " << lastWasSuperior << endl;
    if (lastWasSuperior) {
        position = ceil((start + end) / 2.0);
    } else {
        position = (start + end) / 2 + 1;
    }
//    position = (start + end) / 2;
    return false;


}

/**
 * Add a kmer to a bucket at the right place
 * @param kmer the kmer to add to the bucket
 */
void Bucket::addKmer(Kmer kmer) {
    int position;
    if (find(kmer, position)) { // already in the bucket
        return;
    } else {
        auto itPos = orderedList.begin() + position;
        // We build a superKmer from the kmer and put it in the bucket
        Minimiser minimiserKmer = Minimiser(alpha, this->minimiserLength, kmer);
        Kmer withoutMinimiser = kmer.removePart(minimiserKmer.getPos(), minimiserLength);
        interleavedOrder(withoutMinimiser, minimiserKmer.getPos());

        TYPE prefixLen = minimiserKmer.getPos();
        TYPE suffixLen = kmer.getLength() - minimiserKmer.getPos() - this->minimiserLength;
        int maxLen = (prefixLen < suffixLen) ? suffixLen : prefixLen;
        int fixSize = ceil(log2(kmerLength - minimiserLength + 1));
        SuperKmer toInsert({0});
        toInsert.setBits(0, fixSize, prefixLen);
        toInsert.setBits(fixSize, fixSize, suffixLen);
        toInsert.setBits(2 * fixSize, 4 * maxLen, withoutMinimiser.getValue());

        orderedList.insert(itPos, toInsert);
    }
}

/**
 * Build a Kmer from a superKmer
 * @param superKmer the starting superKmer
 * @return the built Kmer
 */
Kmer Bucket::SKtoKmer(SuperKmer superKmer) {
    int fixBitSize = ceil(log2(kmerLength - minimiserLength + 1));
    uint64_t currentPrefixLen = superKmer.accessBits(0, fixBitSize); // superkmer's prefix length
    uint64_t currentSuffixLen = superKmer.accessBits(fixBitSize, 2 * fixBitSize); // superkmer's suffix length
    uint64_t kmerValue = 0;
    for (int i = currentPrefixLen - 1; i >= 0; i--) {
        int readStart = 2 * fixBitSize + i * 4 + 2;
        kmerValue = (kmerValue << 2) + superKmer.accessBits(readStart, readStart + 2);
    }
    kmerValue = (kmerValue << minimiserLength * 2) + minimiser;
    for (uint64_t i = 0; i < currentSuffixLen; i++) {
        int readStart = 2 * fixBitSize + i * 4;
        kmerValue = (kmerValue << 2) + superKmer.accessBits(readStart, readStart + 2);
    }
    return Kmer(kmerValue, currentPrefixLen + minimiserLength + currentSuffixLen);
}

/**
 * Add a superKmer to a bucket conserving the order
 * @param superKmer the superKmer to add to the bucket
 */
void Bucket::addSuperKmer(const SuperKmer& superKmer) {
    if (orderedList.empty()) {
        return addToList(superKmer);
    }
    // We build a kmer from the superKmer and search it to find the correct position then add it
    Kmer toSearchForPos = SKtoKmer(superKmer);
//    cout << "we testing : " << toSearchForPos.toString() << endl;
    int position;
    if (find(toSearchForPos, position)) return; // Already in
    auto itPos = orderedList.begin() + position;
    orderedList.insert(itPos, superKmer);
}

/**
 * Build a bucket union of two buckets
 * @param toAdd the bucket to join
 * @return a new bucket containing every superKmers of the starting bucket in order without no duplicate
 */
 // Pour opti -> tours lineaires multiples
Bucket Bucket::operator|(const Bucket &toAdd) {
    if (toAdd.minimiser != minimiser || toAdd.kmerLength != kmerLength || toAdd.minimiserLength != minimiserLength) {
        throw std::runtime_error("Error: incompatible buckets");
    }
    Bucket result = Bucket(minimiserLength, minimiser, kmerLength);
    if (orderedList.size() < toAdd.orderedList.size()) {
        result.orderedList = std::vector<SuperKmer>(toAdd.orderedList);
        for (auto &it : orderedList) {
            result.addSuperKmer(it);
        }
        return result;
    } else {
        result.orderedList = std::vector<SuperKmer>(orderedList);
        for (auto &it : toAdd.orderedList) {
            result.addSuperKmer(it);
        }
        return result;
    }
}

/**
 * Getter for the size of the bucket
 * @return the number of superKmer in the bucket
 */
uint64_t Bucket::getListSize() {
    return orderedList.size();
}

/**
 * Getter for the superKmer list of the bucket
 * @return the list of the superKmer in the bucket
 */
std::vector<SuperKmer> Bucket::getListCopy() {
    return std::vector<SuperKmer>(orderedList);
}

/**
 * Getter for the minimiser length of the bucket
 * @return minimiserLength
 */
int Bucket::getMinimiserLength() {
    return minimiserLength;
}

/**
 * Getter for the minimiser value of the bucket
 * @return minimiser
 */
uint64_t Bucket::getMinimiser() {
    return minimiser;
}

/**
 * Getter for the kmerLength of the bucket
 * @return kmerLength
 */
int Bucket::getKmerLength() {
    return minimiserLength;
}

/**
 * Print every superKmer of the bucket as a bitset
 */
void Bucket::print() {
    cout << "content : " << endl;
    for (auto &it : orderedList) {
        it.print();
    }
}

Bucket Bucket::operator&(const Bucket &toIntersect) {
    return Bucket(0, 0, 0);
}


