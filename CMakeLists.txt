cmake_minimum_required(VERSION 3.5...3.28)

project(
  rsiscool
  VERSION 0.0.1
  DESCRIPTION "RS Is Cool by Hamilton Wang")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # clangd

set(CMAKE_BUILD_TYPE Release)
set(BUILD_TYPE "library")

# libraries
set(BREW_CELLAR /opt/homebrew/opt)
set(BREW_INCLUDE /opt/homebrew/include)
set(BREW_LIB /opt/homebrew/lib)

set(DOCTEST_PATH ${BREW_CELLAR}/doctest)
set(EMSCRIPTEN_PATH ${BREW_CELLAR}/emscripten)

# app
set(LIB_SRC src/ReedSolomon.cpp)

include(srcs.cmake)

add_executable(rsiscool src/wasm.cpp ${LIB_SRC} ${NTL_SRC} ${GF2X_SRC})
target_include_directories(
  rsiscool
  SYSTEM
  PRIVATE
  deps/ntl/include
  deps/gf2x
  deps/gf2x/fft
  ${EMSCRIPTEN_PATH}/libexec/system/include)

target_compile_options(rsiscool PRIVATE -O2)
set_target_properties(
  rsiscool
  PROPERTIES
    LINK_FLAGS
    "-s ENVIRONMENT=node -s SINGLE_FILE=0 -s MODULARIZE=1 -s EXPORT_ES6=1 -s 'EXPORT_NAME=rsiscool' --bind --closure 1"
)

add_executable(rsiscool-tests src/test.cpp ${LIB_SRC} ${NTL_SRC} ${NTL_SRC_TEST} ${GF2X_SRC})
set(doctest_DIR ${DOCTEST_PATH}/lib/cmake/doctest)

find_package(doctest REQUIRED)
target_link_libraries(rsiscool-tests PRIVATE doctest::doctest)

# manually edit ntl config.h to disable gmp, cba to add all of gmps source files too
target_include_directories(
  rsiscool-tests
  SYSTEM
  PRIVATE
  ${BREW_INCLUDE} ${BREW_INCLUDE}/gf2x ${DOCTEST_PATH}/include)

target_compile_options(
  rsiscool-tests PRIVATE -O0 -g -Wall -Wextra -Wshadow -pedantic
                   -Wno-missing-field-initializers -Wno-deprecated-declarations)
