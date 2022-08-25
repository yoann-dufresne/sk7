#ifndef SK7_BUCKETMAP_HPP
#define SK7_BUCKETMAP_HPP

#include <unordered_map>

#include "Bucket.hpp"

namespace sk7 {

    class BucketMap {

    public:

        /// Attribute
        std::unordered_map<uint64_t, sk7::Bucket_> *map;

        /// Constructor
        BucketMap();

        /// Requests
        void addKmer(Kmer kmer);
        void addBucket(sk7::Bucket_ bucket);
        bool find(Kmer kmer, int &position) const;

        /// Destructor
        ~BucketMap();

    };
}

#endif //SK7_BUCKETMAP_HPP
