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

set(GMP_VERSION 6.3.0)
set(GMP_PATH /opt/homebrew/Cellar/gmp/${GMP_VERSION})
find_library(
  GMP_LIB
  NAMES gmp
  PATHS ${GMP_PATH}/lib)

set(GIVARO_PATH deps/givaro)
find_library(
  GIVARO_LIB
  NAMES givaro
  PATHS ${GIVARO_PATH}/lib)


# app
set(APP_SRC src/main.cpp src/ReedSolomon.cpp)
add_executable(qriscool ${APP_SRC})

target_link_libraries(qriscool PRIVATE galois ${GIVARO_LIB} ${GMP_LIB})
target_include_directories(qriscool SYSTEM PRIVATE ${GIVARO_PATH}/include ${GMP_PATH}/include)


target_compile_options(qriscool PRIVATE
    -Wall -Wextra -Wshadow -pedantic
    -Wno-missing-field-initializers
    -Wno-deprecated-declarations
  )
