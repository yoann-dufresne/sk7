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
    /// Attributes
    int minimiserLength;
    uint64_t minimiser;
    int kmerLength;
    std::vector<SuperKmer> orderedList;

public:
    /// Constructor
    Bucket(uint64_t minimiserValue);

    /// Add or request methods
    void addToList(SuperKmer superKmer);
    void addKmer(Kmer kmer);
    void addSuperKmer(const SuperKmer& superKmer);
    bool find(Kmer kmer, int &position);

    /// Getter
    uint64_t getListSize();
    std::vector<SuperKmer> getListCopy();

    /// Operator overloading
    Bucket operator|(const Bucket &toAdd); // Union todo : with new compare
    Bucket operator&(Bucket &toIntersect); // Intersection
    Bucket operator^(const Bucket &toXor); // Symmetrical difference todo : union - intersection

    /// Intern functions
    Kmer SKtoKmer(SuperKmer superKmer);
    uint64_t findNextOkPosition(SuperKmer superKmer, std::vector<SuperKmer> list, uint64_t startingPosition);


    /// Misc.
    bool isSorted();
    void print();

};
#endif //SK7_BUCKET_HPP
