# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
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
        FastaReader.cpp
        Kmer.cpp
        Minimiser.cpp
        exampleHash.cpp
        main.cpp
        )

set(HEADERS
        FastaReader.hpp
        Minimiser.hpp
        exampleHash.hpp
        Kmer.hpp
        )

# add the executable
add_executable(sk7 ${SRCS} ${HEADERS})