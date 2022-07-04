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
    if (orderedList.empty()) {
        position = 0;
        return false;
    }
//    cout << endl << endl << "------------------------------" << endl;
//    cout << "searching : " << kmer.toString()  << endl;
    ///PREPARATION OF THE SEARCH
    Minimiser kmerMinimiser = Minimiser(alpha, sk7::m, kmer);
    Kmer withoutMinimiser = kmer.removePart(kmerMinimiser.getPos(), sk7::m);

    uint64_t kmerMask = interleavedOrder(withoutMinimiser, kmerMinimiser.getPos());
//    cout << "striped and interleaved: " << withoutMinimiser.toString() << " of value : " << withoutMinimiser.getValue() << endl;
//    cout << "mask = " << kmerMask << endl;

    int prefixLen = kmerMinimiser.getPos();
    int suffixLen = kmer.getLength() - kmerMinimiser.getPos() - sk7::m;
    int maxLen;


    int needed; // number of bits to compare
    if (prefixLen >= suffixLen) {
        maxLen = prefixLen;
        needed = maxLen * 4;
    } else {
        maxLen = suffixLen;
        needed = maxLen * 4 - 2;
    }


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

        int currentPrefixLen = current.getPrefixLen(); //Current superkmer's prefix length
        int currentSuffixLen = current.getSuffixLen(); //Current superkmer's suffix length

        int currentMaxLen = max(currentPrefixLen, currentSuffixLen);

        uint64_t currentValue = current.getValue(); // superkmer's value
//        cout << "current sk : "; current.print();
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
                toCompare = ((withoutMinimiser.getValue() >> 2) & knownInfo) /*<< 2*/;
                maskedSK >>= 2;
            }
            else toCompare = withoutMinimiser.getValue() & knownInfo;
        }
        else if (found > needed) {
            if (prefixLen < suffixLen) toCompare = withoutMinimiser.getValue() & (knownInfo >> (found - needed - 2));
            else toCompare = withoutMinimiser.getValue() & (knownInfo >> ((found - needed)));
//            /*maskedSK >>= (found - needed);*/
        } else { // (found < needed)
//            cout << "init value : " << withoutMinimiser.getValue() << " shifted : " << (((withoutMinimiser.getValue() >> (needed - found)) & knownInfo)) << endl;
            toCompare = ((withoutMinimiser.getValue() >> (needed - found)) & knownInfo);
//            if (currentPrefixLen < currentSuffixLen) maskedSK >>=2;

        }

//        cout << "before toCompare : " << toCompare << " maskedSK : " << maskedSK << endl;

        if (currentPrefixLen < prefixLen || currentSuffixLen < suffixLen) {
            goto noInfo;
        }

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
                noInfo:
//                cout << "no info" << endl << endl;
                if (start < end) {
                    start += 1;
                    continue;
                }
                if (start == end) { // End of search via linear search
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
    return false;


}

/**
 * Add a kmer to a bucket at the right place
 * @param kmer the kmer to add to the bucket
 */
void Bucket::addKmer(Kmer kmer) {
    if (kmer.length != sk7::k) {
        return;
    }
    int position;
    if (find(kmer, position)) { // already in the bucket
        return;
    } else {
//        cout << "position find = " << position << endl;
        auto itPos = orderedList.begin() + position;
        // We build a superKmer from the kmer and put it in the bucket
        Minimiser minimiserKmer = Minimiser(alpha, sk7::m, kmer);
        Kmer withoutMinimiser = kmer.removePart(minimiserKmer.getPos(), sk7::m);
        interleavedOrder(withoutMinimiser, minimiserKmer.getPos());

        TYPE prefixLen = minimiserKmer.getPos();
        TYPE suffixLen = kmer.getLength() - minimiserKmer.getPos() - sk7::m;
        int maxLen = max(suffixLen, prefixLen);
        SuperKmer toInsert = SuperKmer();
        toInsert.setBits(0, sk7::fixBitSize, prefixLen);
        toInsert.setBits(sk7::fixBitSize, sk7::fixBitSize, suffixLen);
        toInsert.setBits(2 * sk7::fixBitSize, 4 * maxLen, withoutMinimiser.getValue());

        orderedList.insert(itPos, toInsert);
    }
}

/**
 * Build a Kmer from a superKmer
 * @param superKmer the starting superKmer
 * @return the built Kmer
 */
Kmer Bucket::SKtoKmer(SuperKmer superKmer) {
    int PrefixLen = superKmer.getPrefixLen(); // superkmer's prefix length
    int SuffixLen = superKmer.getSuffixLen(); // superkmer's suffix length
    uint64_t kmerValue = 0;
    for (int i = PrefixLen - 1; i >= 0; i--) {
        int readStart = 2 * sk7::fixBitSize + i * 4 + 2;
        kmerValue = (kmerValue << 2) + superKmer.accessBits(readStart, readStart + 2);
    }
    kmerValue = (kmerValue << minimiserLength * 2) + minimiser;
    for (int i = 0; i < SuffixLen; i++) {
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
//    cout << "SK is : \t"; superKmer.print();
    if (orderedList.empty()) {
        return addToList(superKmer);
    }
    // We add each Kmer in the SuperKmer separately for now
    for (auto &kmer : superKmer.split()) {
//        cout << "We adding : " << SKtoKmer(kmer).toString() << endl;
        addKmer(SKtoKmer(kmer));
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
 * For a given SuperKmer, given list and starting position, linearly find the next position where it could be inserted in the list
 * @param superKmer the SuperKmer to position with only one Kmer in it
 * @param list the list to insert the SuperKmer in
 * @param startingPosition the initial position in the list
 * @return the index of a certain position in the list
 */
uint64_t Bucket::findNextOkPosition(SuperKmer superKmer, std::vector<SuperKmer> list, uint64_t startingPosition) {
    uint64_t current = startingPosition;
    uint64_t res = current;
    for (; current < list.size(); current++) {
//        orderedList.at(current).print();
        std::vector<SuperKmer::logic> comparison = SuperKmer::compareSK(superKmer, orderedList.at(current));
        for (auto &logicalValue : comparison) {
            switch (logicalValue) {
                case SuperKmer::EQUAL:
                case SuperKmer::INFERIOR:
                    return current;
                case SuperKmer::SUPERIOR:
                    res = current;
                    break;
                case SuperKmer::INCOMPARABLE:
                    break;
            }
        }
    }
    return res + 1;
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
//    uint64_t i = 0; // loop index for the current Bucket
//    uint64_t j = 0; // loop index for toAdd
//    while (i < orderedList.size() && j < toAdd.orderedList.size()) {
//        switch (SuperKmer::compareSKPerNucleotides(orderedList.at(i), toAdd.orderedList.at(j))) {
//            case SuperKmer::SUPERIOR:
//                result.addToList(toAdd.orderedList.at(j));
//                j++;
//                break;
//            case SuperKmer::INFERIOR:
//                result.addToList(orderedList.at(i));
//                i++;
//                break;
//            case SuperKmer::EQUAL:
//            case SuperKmer::INCOMPARABLE:
//                // We find the one at its place
//            {
//                uint64_t firstI = findNextOkPosition(orderedList.at(i), toAdd.orderedList, j);
//                uint64_t firstJ = findNextOkPosition(toAdd.orderedList.at(j), orderedList, i);
//                if (firstI < firstJ) {
//                    result.addToList(orderedList.at(i));
//                    i++;
//                } else {
//                    result.addToList(toAdd.orderedList.at(j));
//                    j++;
//                }
//            }
//                break;
//
//        }
//    }
//    if (i == orderedList.size()) { // End of this, add the rest of toAdd
//        for (; j < toAdd.orderedList.size(); j++ ) {
//            result.addToList(toAdd.orderedList.at(j));
//        }
//    } else { // End of toXor, add the rest of this
//        for (; i < orderedList.size(); i++ ) {
//            result.addToList(orderedList.at(i));
//        }
//    }
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
    SuperKmer inter;
    while (i < orderedList.size()) {
        SuperKmer current = orderedList.at(i);
        inter = SuperKmer();
        for (uint64_t j = 0; j < toIntersect.getListSize(); j++) {
            inter = inter | (current & toIntersect.orderedList.at(j));
        }
        if (inter != SuperKmer()) {
            result.addToList(inter);
        }
        i++;
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
//    uint64_t i = 0; // loop index for the current Bucket
//    uint64_t j = 0; // loop index for toXor
//    while (i < orderedList.size() && j < toXor.orderedList.size()) {
//        switch (SuperKmer::compareSKPerNucleotides(orderedList.at(i), toXor.orderedList.at(j))) {
//            case SuperKmer::SUPERIOR :
//                result.addToList(toXor.orderedList.at(j));
//                j++;
//                break;
//            case SuperKmer::INFERIOR :
//                result.addToList(orderedList.at(i));
//                i++;
//                break;
//            case SuperKmer::EQUAL :
//                i++;
//                j++;
//                break;
//            case SuperKmer::INCOMPARABLE:
//                result.addToList(toXor.orderedList.at(j));
//                j++;
//                break;
//
//        }
//    }
//    if (i == orderedList.size()) { // End of this, add the rest of toXor
//        for (; j < toXor.orderedList.size(); j++) {
//            result.addToList(toXor.orderedList.at(j));
//        }
//    } else { // End of toXor, add the rest of this
//        for (; i < orderedList.size(); i++) {
//            result.addToList(orderedList.at(i));
//        }
//    }
    return result;
}

/**
 * Check if the bucket is sorted
 * @return true if the bucket is sorted else false
 */
bool Bucket::isSorted() {
    for (uint64_t i = 0; i < orderedList.size() - 1; i++) {
        std::vector<SuperKmer::logic> comparison = SuperKmer::compareSK(orderedList.at(i), orderedList.at(i + 1));
        for (auto &logicalValue : comparison) {
            switch (logicalValue) {
                case SuperKmer::SUPERIOR:
                case SuperKmer::EQUAL:
                    cout << "not sorted in position : " << i << ", " << i + 1 << " with : " << endl;
                    orderedList.at(i).print();
                    orderedList.at(i + 1).print();
                    return false;
                case SuperKmer::INFERIOR:
                    break;
                case SuperKmer::INCOMPARABLE:
                    for (uint64_t k = i + 2; k < orderedList.size(); k++) {
                        std::vector<SuperKmer::logic> linearSearch = SuperKmer::compareSK(orderedList.at(i), orderedList.at(k));
                        for (auto &logicalValueBis : linearSearch) {
                            switch (logicalValueBis) {
                                case SuperKmer::SUPERIOR:
                                case SuperKmer::EQUAL:
                                    cout << "not sorted in position : " << i << ", " << k << " with : " << endl;
                                    orderedList.at(i).print();
                                    orderedList.at(k).print();
                                    return false;
                                case SuperKmer::INFERIOR:
                                case SuperKmer::INCOMPARABLE:
                                    break;
                            }
                        }
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
