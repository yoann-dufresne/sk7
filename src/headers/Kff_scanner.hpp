#ifndef SK7_KFF_SCANNER_HPP
#define SK7_KFF_SCANNER_HPP

#include "kff_io.hpp"
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
    bool sorted;
    uint64_t k;
    uint64_t m;
    uint64_t max;
    uint64_t data_size;

public:

    /// Constructor
    Kff_scanner(char* path, uint8_t m, bool bucketed = false, bool sorted = false);

    /// Methods
    void preparation();
    sk7::Bucket_ readMinimiserSection();
    BucketMap* readAll();

    /// Destructor
    ~Kff_scanner();
};

#endif //SK7_KFF_SCANNER_HPP
