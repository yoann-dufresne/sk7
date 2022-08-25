#include "BucketMap.hpp"

/**
 * Default constructor
 */
BucketMap::BucketMap() {
    this->map = new std::unordered_map<uint64_t , sk7::Bucket>();
};

/// Requests

/**
 * Find the minimizer of a Kmer and interleave it
 * @param kmer the Kmer to prepare for a request
 * @return the value of the minimizer
 */
uint64_t minimizer_finder(Kmer kmer) {
    Minimiser kmerMinimiser = Minimiser(alpha, kmer);
    return kmerMinimiser.getValue();
}



/**
 * Add a Kmer to the proper Bucket
 * @param kmer the Kmer to add
 */
void BucketMap::addKmer(Kmer kmer) {
    uint64_t minimizer = minimizer_finder(kmer);
    try {
        this->map->at(minimizer).addKmer(kmer);
    } catch (const std::out_of_range &e) {
        sk7::Bucket bucket(minimizer);
        bucket.addKmer(kmer);
        this->addBucket(bucket);
    }
}

/**
 * Add a Bucket to this
 * @param bucket the Bucket to add
 */
void BucketMap::addBucket(sk7::Bucket bucket) {
    this->map->insert(std::pair<uint64_t, sk7::Bucket> (bucket.minimiser, bucket));
}

/**
 * Find a Kmer in the correct bucket
 * @param kmer the Kmer to search
 * @param position a reference to an int to store the found position of kmer in its Bucket
 * @return true if kmer was present, false otherwise
 */
bool BucketMap::find(Kmer kmer, int &position) const {
    try {
        std::cout << "minimizer " << std::bitset<20>(minimizer_finder(kmer)) << std::endl;
        this->map->at(minimizer_finder(kmer)).print();
        return this->map->at(minimizer_finder(kmer)).find(kmer, position);
    } catch (const std::out_of_range &e) {
        std::cout << "minimizer not found : " << Kmer(minimizer_finder(kmer), sk7::m).toString() << std::endl;
        return false;
    }
}

/**
 * Destructor
 */
BucketMap::~BucketMap() {
    map->erase(map->begin(), map->cend());
    delete map;
}

