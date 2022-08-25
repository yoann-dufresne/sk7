#include <iostream>
#include <bitset>

using namespace std;

#include "exampleHash.hpp"
#include "Minimiser.hpp"

/**
 * An hash function based on the order A < C < T < G
 * @param kmer the kmer to hash
 * @param length the length of the wanted minimiser
 * @return the value of the lowest length-mer in kmer
 */
hashPos alpha(Kmer kmer, ushort length) {
    uint128_t whole = kmer.getValue();
    uint64_t result = whole;
    short pos = 0;

    uint64_t mask = (1 << (length * 2)) - 1;
    int nb_bit = (kmer.getLength() - length + 1) * 2; // nb of bit comparable

    for (ushort i = 0; i < nb_bit; i += 2) {
        uint64_t current = (whole >> i) & mask;
        if(result > current) { // new minimizer found by value
            pos = nb_bit - i - 1;
            result = current;
        } else if (result == current) { // two possible choice by value, deciding by position
            int middle = (kmer.getLength() - length);
            uint64_t dist_to_middle = abs(pos - middle);
            short current_pos = nb_bit - i - 1;
            uint64_t current_dist_to_middle = abs(current_pos - middle);
            if (dist_to_middle >= current_dist_to_middle) { // new minimizer by position
                pos = current_pos;
            }
        }
    }
    return {result, static_cast<short>(pos / 2)};
}