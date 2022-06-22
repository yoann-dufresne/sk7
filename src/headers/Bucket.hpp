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

    /// Intern function
    uint64_t buildSKMask(const int &prefixLen, const int &suffixLen);

public:
    /// Constructor
    Bucket(int minimiserLength, uint64_t minimiserValue, int kmerLength);

    /// Add or request methods
    void addToList(SuperKmer superKmer);
    void addKmer(Kmer kmer);
    void addSuperKmer(const SuperKmer& superKmer);
    bool find(Kmer kmer, int &position);

    /// Getter
    uint64_t getListSize();
    std::vector<SuperKmer> getListCopy();

    /// Operator overloading
    Bucket operator|(const Bucket &toAdd); // Union
    Bucket operator&(const Bucket &toIntersect); // Intersection
    Bucket operator^(const Bucket &toXor); // Symmetrical difference


    /// Intern functions
    Kmer SKtoKmer(SuperKmer superKmer);
    enum logic {SUPERIOR, INFERIOR, EQUAL, INCOMPARABLE, ENCOMPASSING, ENCOMPASSED, OVERLAPPING};
    logic compareSK(SuperKmer superKmer1, SuperKmer superKmer2);
    uint64_t findNextOkPosition(const SuperKmer& superKmer, std::vector<SuperKmer> list, uint64_t startingPosition);


    /// Misc.
    bool isSorted();
    void print();

};
#endif //SK7_BUCKET_HPP
