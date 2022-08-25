#ifndef SK7_BUCKET_HPP
#define SK7_BUCKET_HPP

#include <cinttypes>
#include <vector>

#include "SuperKmer.hpp"
#include "Kmer.hpp"
#include "Minimizer.hpp"
#include "exampleHash.hpp"
#include "utils.hpp"

namespace sk7 {

    class Bucket_ {

    private:
        /// Attributes
        int minimiserLength;
        int kmerLength;
        std::vector<SuperKmer> orderedList;
    public:

        /// Attributes
        uint64_t minimiser;

        /// Constructor
        Bucket_(uint64_t minimiserValue);

        Bucket_();

        /// Add or request methods
        void addToList(SuperKmer superKmer);
        void addKmer(Kmer kmer);
        void addSuperKmer(const SuperKmer &superKmer);
        bool find(Kmer kmer, int &position);

        /// Getter
        uint64_t getListSize();

        std::vector<SuperKmer> getListCopy();

        /// Operator overloading
        Bucket_ operator|(const Bucket_ &toAdd); // Union
        Bucket_ operator&(Bucket_ &toIntersect); // Intersection
        Bucket_ operator^(const Bucket_ &toXor); // Symmetrical difference

        /// Intern functions
        Kmer SKtoKmer(SuperKmer superKmer);

        uint64_t findNextOkPosition(SuperKmer superKmer, std::vector<SuperKmer> list, uint64_t startingPosition);

        uint64_t nextKmerIndex(const uint64_t &current, const uint64_t &column) const;

        /// Misc.
        bool isSorted();
        void print();

        /// WIP
        static Bucket_ chainedUnion(Bucket_ bucket1, Bucket_ bucket2);
        static bool compatible(const SuperKmer &SK1, const SuperKmer &SK2);

    };
}
#endif //SK7_BUCKET_HPP
