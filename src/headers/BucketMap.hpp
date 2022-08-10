#ifndef SK7_BUCKETMAP_HPP
#define SK7_BUCKETMAP_HPP

#include <unordered_map>

#include "Bucket.hpp"

class BucketMap {

public:

    /// Attribute
    std::unordered_map<uint64_t , sk7::Bucket>* map;

    /// Constructor
    BucketMap();

    /// Requests
    void addKmer(Kmer kmer);
    void addBucket(sk7::Bucket bucket);
    bool find(Kmer kmer, int &position);

    /// Destructor
    ~BucketMap();

};


#endif //SK7_BUCKETMAP_HPP
