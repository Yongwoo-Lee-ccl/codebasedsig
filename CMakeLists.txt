cmake_minimum_required(VERSION 3.12)

set(CMAKE_BUILD_TYPE Release)

SET (CMAKE_C_FLAGS_INIT                "-Wall -std=c99")
SET (CMAKE_C_FLAGS_DEBUG_INIT          "-g")
SET (CMAKE_C_FLAGS_MINSIZEREL_INIT     "-Os -DNDEBUG")
SET (CMAKE_C_FLAGS_RELEASE_INIT        "-O3 -DNDEBUG")
SET (CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")

SET (CMAKE_CXX_FLAGS_INIT                "-Wall")
SET (CMAKE_CXX_FLAGS_DEBUG_INIT          "-g")
SET (CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "-Os -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELEASE_INIT        "-O3 -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")


project(Modified_pqsigRM_Test)

# add_executable(pqsigRM_Test PQCgenKAT_sign.c)
add_executable(pqsigRM_Test main.c)

add_subdirectory(src)
find_package(OpenSSL)

set_target_properties(FILE PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(pqsigRM_Test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(FILE ${CURRENT_LIST_DIR})
target_link_libraries(pqsigRM_Test PRIVATE FILE m OpenSSL::Crypto)
