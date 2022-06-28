#include "headers/Bucket.hpp"

#include <iostream>
#include <cmath>

using namespace std;

/**
 * Initialise a Bucket
 * @param minimiserValue the value of the minimiser
 */
Bucket::Bucket(uint64_t minimiserValue) {
    this->minimiserLength = sk7::m;
    this->minimiser = minimiserValue;
    this->orderedList = std::vector<SuperKmer>();
    this->kmerLength = sk7::k;
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

    int SKheader = 2 * sk7::fixBitSize;


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

        int currentPrefixLen = current.accessBits(0, sk7::fixBitSize); //Current superkmer's prefix length
        int currentSuffixLen = current.accessBits(sk7::fixBitSize, SKheader); //Current superkmer's suffix length
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
//   //         maskedSK >>= 2; //could be used to remove excess 0b00
        } else {
            found = 4 * i;
        }
//        cout << "known info : " << knownInfo << endl;


        /// We align the values
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
            else toCompare = withoutMinimiser.getValue() & (knownInfo >> ((found - needed)));
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
    if (kmer.length != kmerLength) {
        return;
    }
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
    uint64_t PrefixLen = superKmer.accessBits(0, sk7::fixBitSize); // superkmer's prefix length
    uint64_t SuffixLen = superKmer.accessBits(sk7::fixBitSize, 2 * sk7::fixBitSize); // superkmer's suffix length
    uint64_t kmerValue = 0;
    for (int i = PrefixLen - 1; i >= 0; i--) {
        int readStart = 2 * sk7::fixBitSize + i * 4 + 2;
        kmerValue = (kmerValue << 2) + superKmer.accessBits(readStart, readStart + 2);
    }
    kmerValue = (kmerValue << minimiserLength * 2) + minimiser;
    for (uint64_t i = 0; i < SuffixLen; i++) {
        int readStart = 2 * sk7::fixBitSize + i * 4;
        kmerValue = (kmerValue << 2) + superKmer.accessBits(readStart, readStart + 2);
    }
    return Kmer(kmerValue, PrefixLen + minimiserLength + SuffixLen);
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
    int position;
    if (find(toSearchForPos, position)) return; // Already in
    auto itPos = orderedList.begin() + position;
    orderedList.insert(itPos, superKmer);
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
 * For a given SuperKmer and a given list, linearly find the next position where it could be inserted in the list
 * @param superKmer the SuperKmer to position
 * @param list the list to insert the SuperKmer in
 * @param startingPosition the initial position in the list
 * @return the index of a certain position in the list
 */
uint64_t Bucket::findNextOkPosition(const SuperKmer& superKmer, std::vector<SuperKmer> list, uint64_t startingPosition) {
    uint64_t current = startingPosition;
    for (; current < list.size(); current++) {
        switch (SuperKmer::compareSK(superKmer, list.at(current))) {
            case SuperKmer::EQUAL:
            case SuperKmer::ENCOMPASSING:
            case SuperKmer::ENCOMPASSED:
            case SuperKmer::INFERIOR:
                return current;
            case SuperKmer::SUPERIOR:
            case SuperKmer::INCOMPARABLE:
//            case SuperKmer::OVERLAPPING:
                break;
        }
    }
    return current;
}

/**
 * Build a bucket union of two buckets
 * @param toAdd the bucket to join
 * @return a new bucket containing every superKmers of the starting bucket in order without no duplicate
 */
Bucket Bucket::operator|(const Bucket &toAdd) {
    if (toAdd.minimiser != minimiser || toAdd.kmerLength != kmerLength || toAdd.minimiserLength != minimiserLength) {
        throw std::runtime_error("Error: incompatible buckets");
    }
    Bucket result = Bucket(minimiser);
    uint64_t i = 0; // loop index for the current Bucket
    uint64_t j = 0; // loop index for toAdd
    while (i < orderedList.size() && j < toAdd.orderedList.size()) {
        switch (SuperKmer::compareSK(orderedList.at(i), toAdd.orderedList.at(j))) {
            case SuperKmer::SUPERIOR:
                result.addToList(toAdd.orderedList.at(j));
                j++;
                break;
            case SuperKmer::INFERIOR:
                result.addToList(orderedList.at(i));
                i++;
                break;
            case SuperKmer::EQUAL:
            case SuperKmer::ENCOMPASSING:
                result.addToList(orderedList.at(i));
                i++;
                j++;
                break;
            case SuperKmer::ENCOMPASSED:
                result.addToList(orderedList.at(j));
                i++;
                j++;
                break;
            case SuperKmer::INCOMPARABLE:
                // We find the one at its place
            {
                uint64_t firstI = findNextOkPosition(orderedList.at(i), toAdd.orderedList, j);
                uint64_t firstJ = findNextOkPosition(toAdd.orderedList.at(j), orderedList, i);
                if (firstI < firstJ) {
                    result.addToList(orderedList.at(i));
                    i++;
                } else {
                    result.addToList(toAdd.orderedList.at(j));
                    j++;
                }
            }
                break;

        }
    }
    if (i == orderedList.size()) { // End of this, add the rest of toAdd
        for (; j < toAdd.orderedList.size(); j++ ) {
            result.addToList(toAdd.orderedList.at(j));
        }
    } else { // End of toXor, add the rest of this
        for (; i < orderedList.size(); i++ ) {
            result.addToList(orderedList.at(i));
        }
    }
    return result;
}

/**
 * Build the intersection between two Buckets
 * @param toIntersect the bucket to intersect
 * @return the intersection of the two Buckets
 */
Bucket Bucket::operator&(Bucket &toIntersect) {
    if (toIntersect.minimiser != minimiser || toIntersect.kmerLength != kmerLength || toIntersect.minimiserLength != minimiserLength) {
        throw std::runtime_error("Error: incompatible buckets");
    }
    Bucket result = Bucket(minimiser);
    uint64_t i = 0; // loop index for the current Bucket
    uint64_t j = 0; // loop index for toIntersect
    SuperKmer inter;
    while (i < orderedList.size() && j < toIntersect.orderedList.size()) {
        cout << "\tOn compare : " << endl;
        cout << "\t\t"; orderedList.at(i).print();
        cout << "\t\t"; toIntersect.orderedList.at(j).print();
        cout << "\t\t" << SuperKmer::compareSK(orderedList.at(i), toIntersect.orderedList.at(j)) << endl;
        switch (SuperKmer::compareSK(orderedList.at(i), toIntersect.orderedList.at(j))) {
            case SuperKmer::SUPERIOR :
                inter = orderedList.at(i) & toIntersect.orderedList.at(j);
                if (inter != SuperKmer()) {
                    result.addToList(inter);
                }
                j++;
                break;
            case SuperKmer::INFERIOR :
                inter = orderedList.at(i) & toIntersect.orderedList.at(j);
                if (inter != SuperKmer()) {
                    result.addToList(inter);
                }
                i++;
                break;
            case SuperKmer::EQUAL :
                result.addToList(orderedList.at(i));
                i++;
                j++;
                break;
            case SuperKmer::INCOMPARABLE:
                j++;
                break;
            case SuperKmer::ENCOMPASSING:
                result.addToList(toIntersect.orderedList.at(j));
                j++;
                i++;
                break;
            case SuperKmer::ENCOMPASSED:
                result.addToList(orderedList.at(i));
                i++;
                j++;
                break;
//            case SuperKmer::OVERLAPPING:
//                result.addToList(orderedList.at(i) & toIntersect.orderedList.at(j));
//                i++;
//                j++;
//                break;
        }
    }

    return result;
}



/**
 * Build the symmetrical difference between two Buckets
 * @param toXor the Bucket to wor with
 * @return the symmetrical difference between the two Buckets
 */
Bucket Bucket::operator^(const Bucket &toXor) {
    if (toXor.minimiser != minimiser || toXor.kmerLength != kmerLength || toXor.minimiserLength != minimiserLength) {
        throw std::runtime_error("Error: incompatible buckets");
    }
    Bucket result = Bucket(minimiser);
    uint64_t i = 0; // loop index for the current Bucket
    uint64_t j = 0; // loop index for toXor
    while (i < orderedList.size() && j < toXor.orderedList.size()) {
        switch (SuperKmer::compareSK(orderedList.at(i), toXor.orderedList.at(j))) {
            case SuperKmer::SUPERIOR :
                result.addToList(toXor.orderedList.at(j));
                j++;
                break;
            case SuperKmer::INFERIOR :
                result.addToList(orderedList.at(i));
                i++;
                break;
            case SuperKmer::EQUAL :
                i++;
                j++;
                break;
            case SuperKmer::INCOMPARABLE:
                result.addToList(toXor.orderedList.at(j));
                j++;
                break;
            case SuperKmer::ENCOMPASSING:
                break;
            case SuperKmer::ENCOMPASSED:
                break;
//            case SuperKmer::OVERLAPPING:
//                break;
        }
    }
    if (i == orderedList.size()) { // End of this, add the rest of toXor
        for (; j < toXor.orderedList.size(); j++) {
            result.addToList(toXor.orderedList.at(j));
        }
    } else { // End of toXor, add the rest of this
        for (; i < orderedList.size(); i++) {
            result.addToList(orderedList.at(i));
        }
    }
    return result;
}

/**
 * Check if the bucket is sorted
 * @return true if the bucket is sorted else false
 */
bool Bucket::isSorted() {
    for (uint64_t i = 0; i < orderedList.size() - 1; i++) {
        switch (SuperKmer::compareSK(orderedList.at(i), orderedList.at(i + 1))) {
            case SuperKmer::EQUAL:
            case SuperKmer::SUPERIOR:
            case SuperKmer::ENCOMPASSING:
            case SuperKmer::ENCOMPASSED:
                return false;
            case SuperKmer::INFERIOR:
                break;
            case SuperKmer::INCOMPARABLE:
                uint64_t j = i;
                while (j < orderedList.size() - 1) {
                    switch (SuperKmer::compareSK(orderedList.at(j), orderedList.at(j + 1))) {
                        case SuperKmer::SUPERIOR:
                        case SuperKmer::EQUAL:
                        case SuperKmer::ENCOMPASSING:
                        case SuperKmer::ENCOMPASSED:
                            return false;
                        case SuperKmer::INFERIOR:
                            break;
                        case SuperKmer::INCOMPARABLE:
                            j++;
                            continue;
                    }
                    break;
                }
        }
    }
    return true;
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
