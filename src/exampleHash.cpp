#include <bitset>
#include <cmath>

#include "headers/exampleHash.hpp"

/**
 * An hash function based on the alphabetical order
 * @param kmer the kmer to hash
 * @param size the size of the wanted minimiser
 * @return the value of the lowest size-mer in kmer
 */
uint64_t alpha(Kmer kmer, ushort size) {
    uint64_t whole = kmer.getValue();
    uint64_t result = whole;

    uint64_t mask = (1 << (size * 2)) - 1;

    for (ushort i = 0; i < (kmer.getSize() - size + 1) * 2; i += 2) {
        uint64_t current = (whole >> i) & mask;
        result = (result < current) ? result : current;
    }
    return result;
}