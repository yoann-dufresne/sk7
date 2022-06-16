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
    void setBits(const int &start, const int &length, const uint64_t &value);
    void print();
};

#endif //SK7_SUPERKMER_HPP
