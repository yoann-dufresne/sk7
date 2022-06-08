#include <iostream>

#include "exampleHash.hpp"
#include "Minimiser.hpp"

/**
 * An hash function based on the order A < C < T < G
 * @param kmer the kmer to hash
 * @param length the length of the wanted minimiser
 * @return the value of the lowest length-mer in kmer
 */
hashPos alpha(Kmer kmer, ushort length) {
    uint64_t whole = kmer.getValue();
    uint64_t result = whole;
    short pos = 0;

    uint64_t mask = (1 << (length * 2)) - 1;
    int nb_bit = (kmer.getLength() - length + 1) * 2;

    for (ushort i = 0; i < nb_bit; i += 2) {
        uint64_t current = (whole >> i) & mask;
        if(result > current) {
            result = current;
            pos = nb_bit - i - 1;
        }
    }
    return {result, static_cast<short>(pos / 2)};
}