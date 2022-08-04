#ifndef SK7_BUCKET_HPP
#define SK7_BUCKET_HPP

#include <cinttypes>
#include <vector>

#include "SuperKmer.hpp"
#include "Kmer.hpp"
#include "Minimiser.hpp"
#include "exampleHash.hpp"
#include "utils.hpp"

namespace sk7 {

    class Bucket {

    private:
        /// Attributes
        int minimiserLength;
        int kmerLength;
        std::vector<SuperKmer> orderedList;
    public:

        /// Attributes
        uint64_t minimiser;

        /// Constructor
        Bucket(uint64_t minimiserValue);

        Bucket();

        /// Add or request methods
        void addToList(SuperKmer superKmer);

        void addKmer(Kmer kmer);

        void addSuperKmer(const SuperKmer &superKmer);

        bool find(Kmer kmer, int &position);

        /// Getter
        uint64_t getListSize();

        std::vector<SuperKmer> getListCopy();

        /// Operator overloading
        Bucket operator|(const Bucket &toAdd); // Union
        Bucket operator&(Bucket &toIntersect); // Intersection
        Bucket operator^(const Bucket &toXor); // Symmetrical difference

        /// Intern functions
        Kmer SKtoKmer(SuperKmer superKmer);

        uint64_t findNextOkPosition(SuperKmer superKmer, std::vector<SuperKmer> list, uint64_t startingPosition);

        uint64_t nextKmerIndex(const uint64_t &current, const uint64_t &column) const;

        /// Misc.
        bool isSorted();

        void print();

        /// WIP
        static Bucket chainedUnion(Bucket bucket1, Bucket bucket2);

        static bool compatible(const SuperKmer &SK1, const SuperKmer &SK2);

    };
}
#endif //SK7_BUCKET_HPP
