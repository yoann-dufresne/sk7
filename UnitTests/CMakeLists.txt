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

add_executable(unitTest15-5
        unitTest15-5.cpp)

target_link_libraries(unitTest15-5
        sk7)

target_include_directories(unitTest15-5
        PUBLIC
        .
        "${PROJECT_SOURCE_DIR}/extern/lest/include/lest"
        "${PROJECT_SOURCE_DIR}/src/headers"
        )

add_executable(unitTestAltComparison
        unitTestAltComparison.cpp)

target_link_libraries(unitTestAltComparison
        sk7)

target_include_directories(unitTestAltComparison
        PUBLIC
        .
        "${PROJECT_SOURCE_DIR}/extern/lest/include/lest"
        "${PROJECT_SOURCE_DIR}/src/headers"
        )
