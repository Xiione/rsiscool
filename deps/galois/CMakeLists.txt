cmake_minimum_required(VERSION 3.5...3.28)

project(galois
        VERSION 0.0.1
        DESCRIPTION "Galois Field Arithmetic Library by Arash Partow")

add_library(galois STATIC
    src/galois/GaloisField.cpp
    src/galois/GaloisFieldElement.cpp
    src/galois/GaloisFieldPolynomial.cpp
)

target_link_libraries(galois PRIVATE eigen)

target_include_directories(
    galois PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)
