#ifndef SK7_BUCKET_HPP
#define SK7_BUCKET_HPP

#include <cinttypes>
#include <vector>

#include "SuperKmer.hpp"
#include "Kmer.hpp"
#include "Minimiser.hpp"
#include "exampleHash.hpp"
#include "utils.hpp"

class Bucket {

private:
    int minimiserLength;
    uint64_t minimiser;
    int kmerLength;
    std::vector<SuperKmer> orderedList;

public:

    Bucket(int minimiserLength, uint64_t minimiserValue, int kmerLength);
    void addToList(SuperKmer superKmer);
    int isIn(Kmer kmer);
};
#endif //SK7_BUCKET_HPP
