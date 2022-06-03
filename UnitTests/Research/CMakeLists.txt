# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Debug)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

include_directories(${CMAKE_BINARY_DIR})

set(SRCS
        ../../src/SuperKmer.cpp
        ../../src/Kmer.cpp
        ../../src/Bucket.cpp
        ../../src/Minimiser.cpp
        ../../src/utils.cpp
        ../../src/exampleHash.cpp
        testResearch.cpp
        )

set(HEADERS
        ../../src/headers/SuperKmer.hpp
        ../../src/headers/Bucket.hpp
        ../../src/headers/Minimiser.hpp
        ../../src/headers/Kmer.hpp
        ../../src/headers/utils.hpp
        ../../src/headers/exampleHash.hpp
        )

# add the executable
add_executable(testResearch ${SRCS} ${HEADERS})