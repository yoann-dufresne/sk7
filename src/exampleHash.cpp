#include <iostream>

#include "headers/exampleHash.hpp"
#include "headers/Minimiser.hpp"

/**
 * An hash function based on the alphabetical order
 * @param kmer the kmer to hash
 * @param size the size of the wanted minimiser
 * @return the value of the lowest size-mer in kmer
 */
hashPos alpha(Kmer kmer, ushort size) {
    uint64_t whole = kmer.getValue();
    uint64_t result = whole;
    short pos = 0;

    uint64_t mask = (1 << (size * 2)) - 1;
    int nb_bit = (kmer.getSize() - size + 1) * 2;

    for (ushort i = 0; i < nb_bit; i += 2) {
        uint64_t current = (whole >> i) & mask;
        if(result > current) {
            result = current;
            pos = nb_bit - i - 1;
        }
    }
    return {result, static_cast<short>(pos / 2)};
}