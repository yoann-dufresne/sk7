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
uint64_t minimizer_finder(Kmer& kmer) {
    Minimiser kmerMinimiser = Minimiser(alpha, sk7::m, kmer);
//    Kmer withoutMinimiser = kmer.removePart(kmerMinimiser.getPos(), sk7::m);
//    interleavedOrder(withoutMinimiser, kmerMinimiser.getPos());
    return kmerMinimiser.getValue();
}

/**
 * Add a Kmer to the proper Bucket
 * @param kmer the Kmer to add
 */
void BucketMap::addKmer(Kmer kmer) {
    this->map->at(minimizer_finder(kmer)).addKmer(kmer);
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
bool BucketMap::find(Kmer kmer, int &position) {
    this->map->at(minimizer_finder(kmer)).find(kmer, position);
}

/**
 * Destructor
 */
BucketMap::~BucketMap() {
    delete map;
}

