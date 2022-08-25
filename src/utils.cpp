#include "utils.hpp"


using namespace std;
namespace sk7 {

/**
* Build a Kmer in interleaved order from an existing kmer and its minimiser position
* @param kmer striped from its minimiser via Kmer.removePart, will be modify by the function.
* @param minimiserPos the initial position of the minimiser
* @return the mask to find the interleave of the modify kmer
*/
uint64_t interleavedOrder(Kmer &kmer, int minimiserPos) {
    uint64_t finalValue = 0;
    uint64_t mask = 0;
    int prefixLen = minimiserPos;
    int suffixLen = kmer.getLength() - prefixLen;
    int min;
    int max;
    bool prefixIsSmaller;

    if (suffixLen < prefixLen) {
        min = suffixLen;
        max = prefixLen;
        prefixIsSmaller = false;
    } else {
        min = prefixLen;
        max = suffixLen;
        prefixIsSmaller = true;
    }

    int i = 0;
    for (; i < min; i++) { //We read both in the prefix and in the suffix
        finalValue = (finalValue << 2) + kmer.getSubKmer(minimiserPos + i, minimiserPos + i).getValue();
        finalValue = (finalValue << 2) + kmer.getSubKmer(minimiserPos - i - 1, minimiserPos - i - 1).getValue();
        mask = (mask << 4) + 0b1111;
    }

    if (prefixIsSmaller) { //We recover the end of the suffix
        for (; i < max; i++) {
            finalValue = (finalValue << 2) +
                         kmer.getSubKmer(minimiserPos + i, minimiserPos + i).getValue();
            finalValue <<= 2;
            mask = (mask << 4) + 0b1100;
        }

    } else { //We recover the end of the prefix
        for (; i < max; i++) {
            finalValue <<= 2;
            finalValue = (finalValue << 2)
                         + kmer.getSubKmer(minimiserPos - i - 1, minimiserPos - i - 1).getValue();
            mask = (mask << 4) + 0b0011;
        }
    }

    kmer = Kmer(finalValue, kmer.getLength());
    return mask;

}

/**
* Compare two Kmers in value
* @param kmer1 the first Kmer
* @param kmer2 the second Kmer
* @return kmer1 < kmer2 in alphabetical order
*/
bool infId(const Kmer &kmer1, const Kmer &kmer2) {
    return kmer1.getValue() < kmer2.getValue();
}

/**
* Compare two Kmers in value
* @param kmer1 the first Kmer
* @param kmer2 the second Kmer
* @return kmer1 == kmer2 in alphabetical order
*/
bool equalId(const Kmer &kmer1, const Kmer &kmer2) {
    return kmer1.getValue() == kmer2.getValue();
}

}