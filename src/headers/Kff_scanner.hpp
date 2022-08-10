#ifndef SK7_KFF_SCANNER_HPP
#define SK7_KFF_SCANNER_HPP

#include "kff_io.hpp"
#include "Bucket.hpp" // sk7 one
#include "bucket.hpp" // kfftl one
#include "compact.hpp"
#include "BucketMap.hpp"

#include <cstring>
#include <filesystem>


class Kff_scanner {

private:

    /// Attributes
    std::string file_path;
    Kff_file* file;
    bool bucketed;
    bool compacted;
    bool sorted;
    uint64_t k;
    uint64_t m;
    uint64_t max;
    uint64_t data_size;

public:

    /// Constructor
    Kff_scanner(char* path, bool bucketed = true, bool compacted = true, bool sorted = true, uint8_t m = 5);

    /// Methods
    void preparation();
    sk7::Bucket readMinimiserSection();
    BucketMap* readAll();

    /// Destructor
    ~Kff_scanner();
};

#endif //SK7_KFF_SCANNER_HPP
