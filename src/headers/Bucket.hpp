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
    uint64_t buildSKMask(const int &prefixLen, const int &suffixLen);
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


    /// Intern functions
    Kmer SKtoKmer(SuperKmer superKmer);
    enum logic {SUPERIOR, INFERIOR, EQUAL, INCOMPARABLE_TOO_MUCH_INFO, INCOMPARABLE_NOT_ENOUGH_INFO};
    logic compareSK(SuperKmer superKmer1, SuperKmer superKmer2);

    /// Misc.
    void print();

};
#endif //SK7_BUCKET_HPP
