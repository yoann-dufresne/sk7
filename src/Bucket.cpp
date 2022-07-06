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
    interleavedOrder(withoutMinimiser, kmerMinimiser.getPos());

//    cout << "striped and interleaved: " << withoutMinimiser.toString() << " of value : " << withoutMinimiser.getValue() << endl;
//    cout << "mask = " << kmerMask << endl;

    int prefixLen = kmerMinimiser.getPos();

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
//        cout << "currentSK = "; current.print();
//        cout << "prefix Len = " << prefixLen << endl;
        Kmer currentKmer = current.readKmer(prefixLen); // current Kmer
//        cout << "current Kmer = " << currentKmer.toString() << endl;

        if (currentKmer.length == 0) { // No corresponding Kmer in the current SuperKmer
//            cout << "no info" << endl << endl;
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

        uint64_t currentValue = currentKmer.getValue();

//        cout << "search value = " << withoutMinimiser.getValue() << endl;
//        cout << "current value = " << currentValue << endl;

        if (withoutMinimiser.getValue() < currentValue) {
//            cout << "before " << endl << endl;
            start = lastStartWithInformation;
            end = lastPositionWithInformation - 1;
            infoFound = true;
            lastWasSuperior = true;
            continue;
        }

        if (withoutMinimiser.getValue() > currentValue) {
//            cout << "after " << endl << endl;
            start = lastPositionWithInformation + 1;
            lastStartWithInformation = start;
            infoFound = true;
            lastWasSuperior = false;
            continue;
        }

        if (withoutMinimiser.getValue() == currentValue) { //Match or not enough information
                position = middle;
                return true;
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
    uint64_t i = 0; // loop index for the current Bucket
    uint64_t j = 0; // loop index for toAdd
    while (i < orderedList.size() && j < toAdd.orderedList.size()) {
        break;
    }
    return result;
}

/**
 * Build the intersection between two Buckets
 * @param toIntersect the bucket to intersect
 * @return the intersection of the two Buckets
 */
Bucket Bucket::operator&(Bucket &toIntersect) {
//    cout << "---------- start --------------" << endl;
    if (toIntersect.minimiser != minimiser || toIntersect.kmerLength != kmerLength || toIntersect.minimiserLength != minimiserLength) {
        throw std::runtime_error("Error: incompatible buckets");
    }
    uint64_t nbColumn = sk7::k - sk7::m + 1;
    Bucket result = Bucket(minimiser); // the final bucket

    std::vector<uint64_t> idx = std::vector<uint64_t>(nbColumn, 0); // The list of index in the column of the first matrix
    std::vector<uint64_t> idxToIntersect = std::vector<uint64_t>(nbColumn, 0); // The list of index in the column of the second matrix

    for(uint64_t j = 0; j < nbColumn; j++) { // find the first Kmer of the column
        idx.at(j) = nextKmerIndex(0, j);
    }

    for(uint64_t j = 0; j < nbColumn; j++) { // find the first Kmer of the column
        idxToIntersect.at(j) = toIntersect.nextKmerIndex(0, j);
    }

    for (uint64_t j = 0; j < nbColumn; j++) {

        int prefixLen = j;
        int suffixLen = sk7::k - sk7::m - j;

        while (idx.at(j) < orderedList.size() && idxToIntersect.at(j) < toIntersect.orderedList.size()) {

//            cout << "\tj = " << j << endl;
//            orderedList.at(idx.at(j)).print();
//            toIntersect.orderedList.at(idxToIntersect.at(j)).print();

            uint64_t currentSKValue = orderedList.at(idx.at(j)).readKmer(prefixLen).getValue();
            uint64_t currentSKValueToIntersect = toIntersect.orderedList.at(idxToIntersect.at(j)).readKmer(prefixLen).getValue();

//            cout << "current SK = "<< currentSKValue << endl;
//            cout << "current SK i = " << currentSKValueToIntersect << endl;

            if (currentSKValue < currentSKValueToIntersect) {
                idx.at(j) = nextKmerIndex(idx.at(j) + 1, j);
            }
            else if (currentSKValue == currentSKValueToIntersect) {
                SuperKmer toAdd = SuperKmer();
                toAdd.setBits(0, sk7::fixBitSize, prefixLen);
                toAdd.setBits(sk7::fixBitSize, sk7::fixBitSize, suffixLen);
                toAdd.setBits(2*sk7::fixBitSize, 4 * max(prefixLen, suffixLen), currentSKValue);
//                cout << "adding : ";
//                toAdd.print();
                result.addToList(toAdd);
                idx.at(j) = nextKmerIndex(idx.at(j) + 1, j);
                idxToIntersect.at(j) = toIntersect.nextKmerIndex(idxToIntersect.at(j) + 1, j);
            }
            else {
                idxToIntersect.at(j) = toIntersect.nextKmerIndex(idxToIntersect.at(j) + 1, j);
            }
        }
    }
//    cout << "-------- end ----------" << endl;
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
 * @return true if the bucket is sorted with no duplicate else false
 */
bool Bucket::isSorted() {
    std::vector<uint64_t> currentKmers = std::vector<uint64_t>(sk7::k - sk7::m + 1, 0);
    for (uint64_t i = 0; i < orderedList.size(); i++) {
        std::vector<SuperKmer> current = orderedList.at(i).split();
        for (int j = 0; j < sk7::k - sk7::m + 1; j++) {
            if (currentKmers.at(j) >= current.at(j).getValue() && current.at(j) != SuperKmer()) {
                return false;
            } else if (current.at(j) != SuperKmer()) {
                currentKmers.at(j) = current.at(j).getValue();
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

/**
 * Find the index of the next non null kmer in a list of split SuperKmer
 * @param current the starting index in the list
 * @param column the column (linked to the minimiser position) to look in
 * @return the index of the next non null kmer in splitList
 */
uint64_t Bucket::nextKmerIndex(const uint64_t &current, const uint64_t &column) {
    uint64_t i = current;
    for(; i < orderedList.size(); i++) {
        if (orderedList.at(i).readKmer(column).length != 0) {
            return i;
        }
    }
    return i;
}
