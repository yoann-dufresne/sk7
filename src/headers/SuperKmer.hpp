#ifndef SK7_SUPERKMER_HPP
#define SK7_SUPERKMER_HPP

#include <cinttypes>
#include <vector>

#define TYPE uint8_t
#define SIZE 8 // For uint8_t

class SuperKmer {
public:
    SuperKmer(std::vector<TYPE> tab);
    std::vector<TYPE> tab;
    uint64_t accessBits(int start, int end);
};

#endif //SK7_SUPERKMER_HPP
