add_compile_options(-Ofast)

add_library(sk7
    FastaReader.cpp
    Kmer.cpp
    exampleHash.cpp
    Minimiser.cpp
    utils.cpp
    SuperKmer.cpp
    Bucket.cpp
    Kff_scanner.cpp)

target_include_directories(sk7
    PUBLIC
    "${PROJECT_SOURCE_DIR}/src/headers/"
)

target_link_libraries(sk7
        kff
        kfftl
)
