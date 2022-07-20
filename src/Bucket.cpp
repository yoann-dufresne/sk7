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

    uint64_t nbColumn = sk7::k - sk7::m + 1;
    Bucket result = Bucket(minimiser); // the final bucket

    std::vector<uint64_t> idx = std::vector<uint64_t>(nbColumn, 0); // The list of index in the column of the first matrix
    std::vector<uint64_t> idxToAdd = std::vector<uint64_t>(nbColumn, 0); // The list of index in the column of the second matrix

    for(uint64_t j = 0; j < nbColumn; j++) { // find the first Kmer of the column
        idx.at(j) = nextKmerIndex(0, j);
        idxToAdd.at(j) = toAdd.nextKmerIndex(0, j);
    }


    for (uint64_t j = 0; j < nbColumn; j++) {

        uint64_t prefixLen = j;
        uint64_t suffixLen = sk7::k - sk7::m - j;

//        SuperKmer lastAdded;

        while (idx.at(j) < orderedList.size() && idxToAdd.at(j) < toAdd.orderedList.size()) {

//            cout << "\tj = " << j << endl;
//            orderedList.at(idx.at(j)).print();
//            toAdd.orderedList.at(idxToAdd.at(j)).print();

            uint64_t currentKmerValue = orderedList.at(idx.at(j)).readKmer(prefixLen).getValue();
            uint64_t currentKmerValueToAdd = toAdd.orderedList.at(idxToAdd.at(j)).readKmer(prefixLen).getValue();

//            cout << "current K = " << currentKmerValue << endl;
//            cout << "current K a = " << currentKmerValueToAdd << endl;

            SuperKmer SK = SuperKmer();
            SK.setBits(0, sk7::fixBitSize, prefixLen);
            SK.setBits(sk7::fixBitSize, sk7::fixBitSize, suffixLen);
            SK.setBits(2*sk7::fixBitSize, 4 * max(prefixLen, suffixLen), min(currentKmerValue, currentKmerValueToAdd));

            result.addToList(SK);

            if (currentKmerValue < currentKmerValueToAdd) {
                idx.at(j) = nextKmerIndex(idx.at(j) + 1, j);
            }
            else if (currentKmerValue == currentKmerValueToAdd) {
                idx.at(j) = nextKmerIndex(idx.at(j) + 1, j);
                idxToAdd.at(j) = toAdd.nextKmerIndex(idxToAdd.at(j) + 1, j);
            }
            else {
                idxToAdd.at(j) = toAdd.nextKmerIndex(idxToAdd.at(j) + 1, j);
            }
        }

        if (idx.at(j) == orderedList.size()) { // end of this' column, add the rest of toAdd's
            while (idxToAdd.at(j) < toAdd.orderedList.size()) {

                uint64_t currentKmerValueToAdd = toAdd.orderedList.at(idxToAdd.at(j)).readKmer(prefixLen).getValue();
                SuperKmer SK = SuperKmer();

                SK.setBits(0, sk7::fixBitSize, prefixLen);
                SK.setBits(sk7::fixBitSize, sk7::fixBitSize, suffixLen);
                SK.setBits(2*sk7::fixBitSize, 4 * max(prefixLen, suffixLen), currentKmerValueToAdd);

                result.addToList(SK);
                idxToAdd.at(j) = toAdd.nextKmerIndex(idxToAdd.at(j) + 1, j);
            }
        } else { // end of toAdd's column, add the rest of this'
            while (idx.at(j) < orderedList.size()) {

                uint64_t currentKmerValue = orderedList.at(idx.at(j)).readKmer(prefixLen).getValue();
                SuperKmer SK = SuperKmer();

                SK.setBits(0, sk7::fixBitSize, prefixLen);
                SK.setBits(sk7::fixBitSize, sk7::fixBitSize, suffixLen);
                SK.setBits(2 * sk7::fixBitSize, 4 * max(prefixLen, suffixLen), currentKmerValue);

                result.addToList(SK);
                idx.at(j) = nextKmerIndex(idx.at(j) + 1, j);
            }
        }
    }
//    cout << "-------- end ----------" << endl;
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
        idxToIntersect.at(j) = toIntersect.nextKmerIndex(0, j);
    }


    for (uint64_t j = 0; j < nbColumn; j++) {

        uint64_t prefixLen = j;
        uint64_t suffixLen = sk7::k - sk7::m - j;

        while (idx.at(j) < orderedList.size() && idxToIntersect.at(j) < toIntersect.orderedList.size()) {

//            cout << "\tj = " << j << endl;
//            orderedList.at(idx.at(j)).print();
//            toIntersect.orderedList.at(idxToIntersect.at(j)).print();

            uint64_t currentKmerValue = orderedList.at(idx.at(j)).readKmer(prefixLen).getValue();
            uint64_t currentKmerValueToIntersect = toIntersect.orderedList.at(idxToIntersect.at(j)).readKmer(prefixLen).getValue();

//            cout << "current SK = "<< currentKmerValue << endl;
//            cout << "current SK i = " << currentKmerValueToIntersect << endl;

            if (currentKmerValue < currentKmerValueToIntersect) {
                idx.at(j) = nextKmerIndex(idx.at(j) + 1, j);
            }
            else if (currentKmerValue == currentKmerValueToIntersect) {
                SuperKmer toAdd = SuperKmer();
                toAdd.setBits(0, sk7::fixBitSize, prefixLen);
                toAdd.setBits(sk7::fixBitSize, sk7::fixBitSize, suffixLen);
                toAdd.setBits(2*sk7::fixBitSize, 4 * max(prefixLen, suffixLen), currentKmerValue);
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
    uint64_t nbColumn = sk7::k - sk7::m + 1;
    Bucket result = Bucket(minimiser); // the final bucket

    std::vector<uint64_t> idx = std::vector<uint64_t>(nbColumn, 0); // The list of index in the column of the first matrix
    std::vector<uint64_t> idxToXor = std::vector<uint64_t>(nbColumn, 0); // The list of index in the column of the second matrix

    for(uint64_t j = 0; j < nbColumn; j++) { // find the first Kmer of the column
        idx.at(j) = nextKmerIndex(0, j);
        idxToXor.at(j) = toXor.nextKmerIndex(0, j);
    }

    for (uint64_t j = 0; j < nbColumn; j++) {

        uint64_t prefixLen = j;
        uint64_t suffixLen = sk7::k - sk7::m - j;

        while (idx.at(j) < orderedList.size() && idxToXor.at(j) < toXor.orderedList.size()) {

//            cout << "\tj = " << j << endl;
//            orderedList.at(idx.at(j)).print();
//            toIntersect.orderedList.at(idxToXor.at(j)).print();

            uint64_t currentKmerValue = orderedList.at(idx.at(j)).readKmer(prefixLen).getValue();
            uint64_t currentKmerValueToXor = toXor.orderedList.at(idxToXor.at(j)).readKmer(prefixLen).getValue();

//            cout << "current SK = "<< currentKmerValue << endl;
//            cout << "current SK i = " << currentKmerValueToXor << endl;

            if (currentKmerValue < currentKmerValueToXor) {
                SuperKmer toAdd = SuperKmer();
                toAdd.setBits(0, sk7::fixBitSize, prefixLen);
                toAdd.setBits(sk7::fixBitSize, sk7::fixBitSize, suffixLen);
                toAdd.setBits(2*sk7::fixBitSize, 4 * max(prefixLen, suffixLen), currentKmerValue);
//                cout << "adding : ";
//                toAdd.print();
                result.addToList(toAdd);
                idx.at(j) = nextKmerIndex(idx.at(j) + 1, j);
            }
            else if (currentKmerValue == currentKmerValueToXor) {
                idx.at(j) = nextKmerIndex(idx.at(j) + 1, j);
                idxToXor.at(j) = toXor.nextKmerIndex(idxToXor.at(j) + 1, j);
                continue;
            }
            else {
                SuperKmer toAdd = SuperKmer();
                toAdd.setBits(0, sk7::fixBitSize, prefixLen);
                toAdd.setBits(sk7::fixBitSize, sk7::fixBitSize, suffixLen);
                toAdd.setBits(2*sk7::fixBitSize, 4 * max(prefixLen, suffixLen), currentKmerValueToXor);
//                cout << "adding : ";
//                toAdd.print();
                result.addToList(toAdd);
                idxToXor.at(j) = toXor.nextKmerIndex(idxToXor.at(j) + 1, j);
            }
        }

        if (idx.at(j) == orderedList.size()) { // end of this' column, add the rest of toAdd's
            while (idxToXor.at(j) < toXor.orderedList.size()) {

                uint64_t currentKmerValueToXor = toXor.orderedList.at(idxToXor.at(j)).readKmer(prefixLen).getValue();
                SuperKmer SK = SuperKmer();

                SK.setBits(0, sk7::fixBitSize, prefixLen);
                SK.setBits(sk7::fixBitSize, sk7::fixBitSize, suffixLen);
                SK.setBits(2*sk7::fixBitSize, 4 * max(prefixLen, suffixLen), currentKmerValueToXor);

                result.addToList(SK);
                idxToXor.at(j) = toXor.nextKmerIndex(idxToXor.at(j) + 1, j);
            }
        } else { // end of toAdd's column, add the rest of this'
            while (idx.at(j) < orderedList.size()) {

                uint64_t currentKmerValue = orderedList.at(idx.at(j)).readKmer(prefixLen).getValue();
                SuperKmer SK = SuperKmer();

                SK.setBits(0, sk7::fixBitSize, prefixLen);
                SK.setBits(sk7::fixBitSize, sk7::fixBitSize, suffixLen);
                SK.setBits(2 * sk7::fixBitSize, 4 * max(prefixLen, suffixLen), currentKmerValue);

                result.addToList(SK);
                idx.at(j) = nextKmerIndex(idx.at(j) + 1, j);
            }
        }
    }
//    cout << "-------- end ----------" << endl;
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
uint64_t Bucket::nextKmerIndex(const uint64_t &current, const uint64_t &column) const {
    uint64_t i = current;
    for(; i < orderedList.size(); i++) {
        if (orderedList.at(i).readKmer(column).length != 0) {
            return i;
        }
    }
    return i;
}


/// TEST ZONE

/**
 * Build the union of two compatible buckets trying to conserve the compression
 * @param bucket1 the first Bucket
 * @param bucket2 the second Bucket
 * @return the built union
 */
Bucket Bucket::chainedUnion(Bucket bucket1, Bucket bucket2) {

    /// Initialisation
    if (bucket1.minimiser != bucket2.minimiser ||
        bucket1.kmerLength != bucket2.kmerLength ||
        bucket1.minimiserLength != bucket2.minimiserLength) {
        throw std::runtime_error("Error: incompatible buckets");
    }

    Bucket result = Bucket(bucket1.minimiser);

    uint64_t nbColumn = sk7::k - sk7::m + 1; // = max number of Kmer per SuperKmer

    std::vector<uint64_t> idx1 = std::vector<uint64_t>(nbColumn, 0); // The list of index in the column of the first matrix
    std::vector<uint64_t> idx2 = std::vector<uint64_t>(nbColumn, 0); // The list of index in the column of the second matrix

    uint64_t current_line = max(bucket1.getListSize(), bucket2.getListSize());
    uint64_t current_column = 0;
    uint64_t current_bucket = 0;

    for(uint64_t column = 0; column < nbColumn; column++) { // find the first Kmer of the column (column number = prefix len)

        uint64_t line1 = bucket1.nextKmerIndex(0, column);
        uint64_t line2 = bucket2.nextKmerIndex(0, column);

        idx1.at(column) = line1;
        idx2.at(column) = line2;

        uint64_t kmer1;
        uint64_t kmer2;

        try {
            kmer1 = bucket1.orderedList.at(line1).extract(column).getValue();
        } catch (...) {
            kmer1 = UINT64_MAX;
        }

        try {
            kmer2 = bucket2.orderedList.at(line2).extract(column).getValue();
        } catch (...) {
            kmer2 = UINT64_MAX;
        }

        if (kmer1 < kmer2 && line1 < current_line) {
            current_line = line1;
            current_column = column;
            current_bucket = 1;
        }

        else if (kmer2 < kmer1 && line2 < current_line) {
            current_line = line2;
            current_column = column;
            current_bucket = 2;
        }
    }

    /// Iterations
    while ((current_bucket == 1)? current_line < bucket1.getListSize() : current_line < bucket2.getListSize()) {

//        cout << "NEW TOUR" << endl;
//        cout << "current line = " << current_line << endl;
//        cout << "current column = " << current_column << endl;
//        cout << "current bucket = " << current_bucket << endl;

        SuperKmer toAdd; // The SuperKmer to add at the end of the iteration
        SuperKmer current; // The chosen Kmer for the iteration
        SuperKmer last; // Last Kmer seen to check for compatibility

        if (current_bucket == 1) { // We take current in the first bucket
            toAdd = current = last = bucket1.orderedList.at(idx1.at(current_column)).extract(current_column);
            idx1.at(current_column) = bucket1.nextKmerIndex(idx1.at(current_column) + 1, current_column);
            try {
                if (bucket2.orderedList.at(idx2.at(current_column)).extract(current_column) == current) { // check for double possibility
                    idx2.at(current_column) = bucket2.nextKmerIndex(idx2.at(current_column) + 1, current_column);
                }
            } catch (...) {

            }

        } else { // We take current in the second bucket
            toAdd = current = last = bucket2.orderedList.at(idx2.at(current_column)).extract(current_column);
            idx2.at(current_column) = bucket2.nextKmerIndex(idx2.at(current_column) + 1, current_column);
            try {
                if (bucket1.orderedList.at(idx1.at(current_column)).extract(current_column) == current) { // check for double possibility
                    idx1.at(current_column) = bucket1.nextKmerIndex(idx1.at(current_column) + 1, current_column);
                }
            } catch (...) {

            }

        }

//        cout << "current = "; current.print();

        for (uint64_t i = current_column + 1; i < nbColumn; i++) { // Looking for possibility in upper columns
//            cout << "going up : " << i << endl;
//            cout << "\tidx1 -> " << idx1.at(i) << " idx2 -> " << idx2.at(i) << endl;
            if (idx1.at(i) < bucket1.getListSize() && idx2.at(i) < bucket2.getListSize()) { // two Kmers to compare
                SuperKmer neighbor1 = bucket1.orderedList.at(idx1.at(i)).extract(i);
                SuperKmer neighbor2 = bucket2.orderedList.at(idx2.at(i)).extract(i);

//                cout << "\tneighbor1 ";
//                neighbor1.print();
//                cout << "\tneighbor2 ";
//                neighbor2.print();

                if (neighbor1.readKmer(i).getValue() < neighbor2.readKmer(i).getValue()
                    && Bucket::compatible(neighbor1, last)) { // Adding the Kmer of the first bucket if compatible
//                    cout << "adding 1" << endl;
                    last = neighbor1;
                    toAdd = toAdd | neighbor1;
                    idx1.at(i) = bucket1.nextKmerIndex(idx1.at(i) + 1, i);
                    continue;

                } else if (neighbor2.readKmer(i).getValue() < neighbor1.readKmer(i).getValue()
                           && Bucket::compatible(neighbor2, last)) { // Adding the Kmer of the second bucket if compatible
//                    cout << "adding 2" << endl;
                    last = neighbor2;
                    toAdd = toAdd | neighbor2;
                    idx2.at(i) = bucket2.nextKmerIndex(idx2.at(i) + 1, i);
                    continue;
                } else if (neighbor2.readKmer(i).getValue() == neighbor1.readKmer(i).getValue()
                           && Bucket::compatible(neighbor2, last)) { // The two Kmers are equals -> adding the two if compatible
//                    cout << "adding 1 & 2" << endl;
                    last = neighbor1;
                    toAdd = toAdd | neighbor2;
                    idx1.at(i) = bucket1.nextKmerIndex(idx1.at(i) + 1, i);
                    idx2.at(i) = bucket2.nextKmerIndex(idx2.at(i) + 1, i);
                    continue;
                }

                break;
            }

            else if (idx2.at(i) < bucket2.getListSize()) { // No more Kmer in this column of bucket1
                SuperKmer neighbor2 = bucket2.orderedList.at(idx2.at(i)).extract(i);
                if (Bucket::compatible(neighbor2, last)) {
                    last = neighbor2;
                    toAdd = toAdd | neighbor2;
                    idx2.at(i) = bucket2.nextKmerIndex(idx2.at(i) + 1, i);
                    continue;
                }
            }

            else if (idx1.at(i) < bucket1.getListSize()) { // No more Kmer in this column of bucket2
                SuperKmer neighbor1 = bucket1.orderedList.at(idx1.at(i)).extract(i);
                if (Bucket::compatible(neighbor1, last)) {
                    last = neighbor1;
                    toAdd = toAdd | neighbor1;
                    idx1.at(i) = bucket1.nextKmerIndex(idx1.at(i) + 1, i);
                    continue;
                }
            }

            else
                break;
        }

        last = current; // reset of last

        for (int64_t i = current_column - 1; i >= 0 ; i--) { // Looking for possibility in lower columns
//            cout << "going down : " << i << endl;
            if (idx1.at(i) < bucket1.getListSize() && idx2.at(i) < bucket2.getListSize()) { // two Kmers to compare
                SuperKmer neighbor1 = bucket1.orderedList.at(idx1.at(i)).extract(i);
                SuperKmer neighbor2 = bucket2.orderedList.at(idx2.at(i)).extract(i);

//                cout << "\tneighbor1 ";
//                neighbor1.print();
//                cout << "\tneighbor2 ";
//                neighbor2.print();

                if (neighbor1.readKmer(i).getValue() < neighbor2.readKmer(i).getValue()
                    && Bucket::compatible(neighbor1, last)) { // Adding the Kmer of the first bucket if compatible
//                    cout << "adding 1" << endl;
                    last = neighbor1;
                    toAdd = toAdd | neighbor1;
                    idx1.at(i) = bucket1.nextKmerIndex(idx1.at(i) + 1, i);
                    continue;
                } else if (neighbor2.readKmer(i).getValue() < neighbor1.readKmer(i).getValue()
                           && Bucket::compatible(neighbor2, last)) { // Adding the Kmer of the second bucket if compatible
//                    cout << "adding 2" << endl;
                    last = neighbor2;
                    toAdd = toAdd | neighbor2;
                    idx2.at(i) = bucket2.nextKmerIndex(idx2.at(i) + 1, i);
                    continue;
                } else if (neighbor2.readKmer(i).getValue() == neighbor1.readKmer(i).getValue()
                           && Bucket::compatible(neighbor2, last)) { // The two Kmers are equals -> adding the two if compatible
//                    cout << "adding 1 & 2" << endl;
                    last = neighbor1;
                    toAdd = toAdd | neighbor2;
                    idx1.at(i) = bucket1.nextKmerIndex(idx1.at(i) + 1, i);
                    idx2.at(i) = bucket2.nextKmerIndex(idx2.at(i) + 1, i);
                    continue;
                }

                break;
            }

            else if (idx2.at(i) < bucket2.getListSize()) { // No more Kmer in this column of bucket1
                SuperKmer neighbor2 = bucket2.orderedList.at(idx2.at(i)).extract(i);
                if (Bucket::compatible(neighbor2, last)) {
                    last = neighbor2;
                    toAdd = toAdd | neighbor2;
                    idx2.at(i) = bucket2.nextKmerIndex(idx2.at(i) + 1, i);
                    continue;
                }
            }

            else if (idx1.at(i) < bucket1.getListSize()) { // No more Kmer in this column of bucket2
                SuperKmer neighbor1 = bucket1.orderedList.at(idx1.at(i)).extract(i);
                if (Bucket::compatible(neighbor1, last)) {
                    last = neighbor1;
                    toAdd = toAdd | neighbor1;
                    idx1.at(i) = bucket1.nextKmerIndex(idx1.at(i) + 1, i);
                    continue;
                }
            }

            else
                break;

        }

//        cout << "\t\t\ttoAdd : " ; toAdd.print();
        result.addToList(toAdd); // Adding the built SuperKmer

        current_line = max(bucket1.getListSize(), bucket2.getListSize());

        for(uint64_t column = 0; column < nbColumn; column++) { // find the coordinates of the next current

            uint64_t line1 = idx1.at(column);
            uint64_t line2 = idx2.at(column);

            uint64_t kmer1;
            uint64_t kmer2;

            try {
                kmer1 = bucket1.orderedList.at(line1).extract(column).getValue();
            } catch (...) {
                kmer1 = UINT64_MAX;
            }

            try {
                kmer2 = bucket2.orderedList.at(line2).extract(column).getValue();
            } catch (...) {
                kmer2 = UINT64_MAX;
            }

            if (kmer1 < kmer2 && line1 < current_line) {
                    current_line = line1;
                    current_column = column;
                    current_bucket = 1;
            }

            else if (kmer2 < kmer1 && line2 < current_line) {
                    current_line = line2;
                    current_column = column;
                    current_bucket = 2;
            }

        }
    }

    return result;
}



bool Bucket::compatible(const SuperKmer& SK1, const SuperKmer& SK2) {
//    cout << endl;
    uint64_t rangeLen = sk7::k - sk7::m - 1;
    uint64_t mask = (1u << rangeLen * 2) - 1;
//    cout << "mask    = " << bitset<64>(mask) << endl;

    uint64_t kmerValue1 = 0;
    uint64_t kmerValue2 = 0;

    if (SK1.getPrefixLen() < SK2.getPrefixLen()) {
        kmerValue1 = (SK1.nonInterleavedKmerValue() >> 2) ;
        kmerValue2 = SK2.nonInterleavedKmerValue();
    }
    else {
        kmerValue1 = SK1.nonInterleavedKmerValue() ;
        kmerValue2 = (SK2.nonInterleavedKmerValue() >> 2);
    }

//    cout << "value 1 = " << bitset<64>(kmerValue1) << endl;
//    cout << "value 2 = " << bitset<64>(kmerValue2) << endl;

    return (kmerValue1 & mask) == (kmerValue2 & mask);
}