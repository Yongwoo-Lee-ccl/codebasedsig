cmake_minimum_required(VERSION 3.12)

set(CMAKE_BUILD_TYPE Debug)

project(Modified_pqsigRM_Test)

# add_executable(pqsigRM_Test PQCgenKAT_sign.c)
add_executable(pqsigRM_Test main.c)

add_subdirectory(src)
find_package(OpenSSL)

set_target_properties(FILE PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(pqsigRM_Test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(FILE ${CURRENT_LIST_DIR})
target_link_libraries(pqsigRM_Test PRIVATE FILE m OpenSSL::Crypto)
