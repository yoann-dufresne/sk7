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
        ../../src/FastaReader.cpp
        ../../src/Kmer.cpp
        ../../src/Minimiser.cpp
        ../../src/exampleHash.cpp
        testKmerMini.cpp
        )

set(HEADERS
        ../../src/headers/FastaReader.hpp
        ../../src/headers/Minimiser.hpp
        ../../src/headers/exampleHash.hpp
        ../../src/headers/Kmer.hpp
        )

add_executable(testKmerMini ${SRCS} ${HEADERS})
