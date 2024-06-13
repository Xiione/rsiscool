cmake_minimum_required(VERSION 3.5...3.28)

project(
  rsiscool
  VERSION 0.0.1
  DESCRIPTION "QR Is Cool by Hamilton Wang")

set(CMAKE_C_COMPILER "emcc")
set(CMAKE_CXX_COMPILER "em++")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # clangd

set(CMAKE_BUILD_TYPE Release)
set(BUILD_TYPE "library")

# libraries
set(BREW_CELLAR /opt/homebrew/opt)
set(BREW_INCLUDE /opt/homebrew/include)

# set(NTL_PATH ${BREW_CELLAR}/ntl)

set(EMSCRIPTEN_PATH ${BREW_CELLAR}/emscripten)

# app
set(LIB_SRC src/ReedSolomon.cpp)

set(NTL deps/ntl/src)
set(NTL_SRC
    ${NTL}/lip.cpp
    ${NTL}/tools.cpp
    ${NTL}/WordVector.cpp
    ${NTL}/ZZ.cpp
    ${NTL}/GF2.cpp
    ${NTL}/GF2X.cpp
    ${NTL}/GF2X1.cpp
    ${NTL}/GF2XVec.cpp
    ${NTL}/GF2E.cpp
    ${NTL}/GF2EX.cpp
    ${NTL}/vec_GF2E.cpp
    ${NTL}/mat_GF2E.cpp)

set(GF2X deps/gf2x)
set(GF2X_SRC
    ${GF2X}/gf2x.c
    ${GF2X}/toom.c
    ${GF2X}/toom-gpl.c
    ${GF2X}/fft/gf2x-ternary-fft.c)

add_executable(rsiscool ${LIB_SRC} ${NTL_SRC} ${GF2X_SRC})
target_include_directories(
  rsiscool SYSTEM PRIVATE deps/ntl/include deps/gf2x deps/gf2x/fft
                          ${EMSCRIPTEN_PATH}/libexec/system/include)

# add_library(NTL_LIB STATIC IMPORTED) set(NTL_PATH /usr/local/lib/libntl.a)
# find_library( NTL_LIB NAMES ntl libntl PATHS /usr/local/lib REQUIRED
# NO_DEFAULT_PATH) set_property(TARGET NTL_LIB PROPERTY IMPORTED_LOCATION
# ${NTL_PATH}) target_link_libraries(rsiscool ${NTL_LIB})

target_compile_options(
  rsiscool
  PRIVATE -O2
          -Wall
          -Wextra
          -Wshadow
          -pedantic
          -Wno-missing-field-initializers
          -Wno-deprecated-declarations)

set_target_properties(
  rsiscool
  PROPERTIES
    LINK_FLAGS
    "-s ENVIRONMENT=node -s SINGLE_FILE=0 -s MODULARIZE=1 -s EXPORT_ES6=1 -s 'EXPORT_NAME=rsiscool' --bind --closure 1"
)
