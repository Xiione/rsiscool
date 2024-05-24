cmake_minimum_required(VERSION 3.5...3.28)

project(
  qriscool
  VERSION 0.0.1
  DESCRIPTION "QR Is Cool by Hamilton Wang")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # clangd

set(APP_SRC src/main.cpp src/ReedSolomon.cpp)

set(GIVARO_DIR deps/givaro)
include_directories(${GIVARO_DIR}/include)
find_library(
  GIVARO_LIB
  NAMES givaro
  PATHS ${GIVARO_LIB}/lib)

add_executable(qriscool ${APP_SRC})
target_link_libraries(qriscool PRIVATE GIVARO_LIB)
