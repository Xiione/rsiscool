cmake_minimum_required(VERSION 3.5...3.28)

project(
  qriscool
  VERSION 0.0.1
  DESCRIPTION "QR Is Cool by Hamilton Wang")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # clangd

# libraries
add_subdirectory(deps/galois)

set(BREW_PATH /opt/homebrew/Cellar)

# set(GMP_VERSION 6.3.0)
# set(GMP_PATH ${BREW_PATH}/gm/${GMP_VERSION})
# find_library(
#   GMP_LIB
#   NAMES gmp
#   PATHS ${GMP_PATH}/lib)
#
# set(GF2X_VERSION 1.3.99)
# set(GF2X_PATH ${BREW_PATH}/gf2x/${GF2X_VERSION})
# find_library(
#   GF2X_LIB
#   NAMES gf2x
#   PATHS ${GF2X_PATH}/lib)

set(NTL_VERSION 11.5.1)
set(NTL_PATH ${BREW_PATH}/ntl/${NTL_VERSION})
find_library(
  NTL_LIB
  NAMES ntl
  PATHS ${NTL_PATH}/lib)

# app
set(APP_SRC src/main.cpp src/ReedSolomon.cpp)
add_executable(qriscool ${APP_SRC})

target_link_libraries(qriscool PRIVATE galois ${NTL_LIB})
target_include_directories(
  qriscool SYSTEM PRIVATE ${NTL_PATH}/include)

target_compile_options(
  qriscool PRIVATE -g -Wall -Wextra -Wshadow -pedantic
                   -Wno-missing-field-initializers -Wno-deprecated-declarations)
