
add_library(sk7
    FastaReader.cpp
    Kmer.cpp
    exampleHash.cpp
    Minimiser.cpp
    utils.cpp
    SuperKmer.cpp
    Bucket.cpp
)

target_include_directories(sk7
    PUBLIC
    "${PROJECT_SOURCE_DIR}/src/headers/"
)
