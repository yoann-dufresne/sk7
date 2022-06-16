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
    /// Constructor
    Bucket(int minimiserLength, uint64_t minimiserValue, int kmerLength);

    /// Add or request methods
    void addToList(SuperKmer superKmer);
    void addKmer(Kmer kmer);
    void addSuperKmer(const SuperKmer& superKmer);
    bool find(Kmer kmer, int &position);

    /// Getter
    int getMinimiserLength();
    uint64_t getMinimiser();
    int getKmerLength();
    uint64_t getListSize();
    std::vector<SuperKmer> getListCopy();

    /// Operator overloading
    Bucket operator|(const Bucket &toAdd); // Union
    Bucket operator&(const Bucket &toIntersect); // Intersection
    Bucket operator^(const Bucket &toXor); // symmetrical difference

    /// Misc.
    Kmer SKtoKmer(SuperKmer superKmer);
    void print();

};
#endif //SK7_BUCKET_HPP
