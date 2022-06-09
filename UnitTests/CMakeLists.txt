add_executable(unitTest
        unitTest.cpp)

target_link_libraries(unitTest
        sk7)

target_include_directories(unitTest
        PUBLIC
        .
        "${PROJECT_SOURCE_DIR}/extern/lest/include/lest"
        "${PROJECT_SOURCE_DIR}/src/headers"
        )