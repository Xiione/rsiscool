cmake_minimum_required(VERSION 3.5...3.28)

project(
  qriscool
  VERSION 0.0.1
  DESCRIPTION "QR Is Cool by Hamilton Wang")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # clangd

# libraries
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

include_directories(${GIVARO_PATH}/include ${GMP_PATH}/include)

# app
set(APP_SRC src/main.cpp src/ReedSolomon.cpp)
add_executable(qriscool ${APP_SRC})
target_link_libraries(qriscool PRIVATE ${GIVARO_LIB} ${GMP_LIB})
