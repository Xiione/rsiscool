cmake_minimum_required(VERSION 3.5...3.28)

project(
  rsiscool
  VERSION 0.0.2
  DESCRIPTION "RS Is Cool by Hamilton Wang")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # clangd

set(CMAKE_BUILD_TYPE Release)
set(BUILD_TYPE "library")

set(BREW_CELLAR /opt/homebrew/opt)
set(EMSCRIPTEN_PATH ${BREW_CELLAR}/emscripten)

# app
set(LIB_SRC src/ReedSolomon.cpp)
add_subdirectory(deps/galois)


add_executable(rsiscool src/wasm.cpp ${LIB_SRC})

include_directories(${EMSCRIPTEN_PATH}/libexec/system/include)
target_link_libraries(rsiscool PRIVATE galois)
target_compile_options(rsiscool PRIVATE -O3 -fno-exceptions)
set_target_properties(
  rsiscool
  PROPERTIES
    LINK_FLAGS "\
    -s ENVIRONMENT=worker\
    -s SINGLE_FILE=0\
    -s MODULARIZE=1\
    -s EXPORT_ES6=1\
    -s 'EXPORT_NAME=rsiscool'\
    --bind --closure 1 -flto\
    "
)

add_executable(rsiscool-tests src/test.cpp ${LIB_SRC})

set(DOCTEST_PATH ${BREW_CELLAR}/doctest)
set(doctest_DIR ${DOCTEST_PATH}/lib/cmake/doctest)
find_package(doctest REQUIRED)
target_link_libraries(rsiscool-tests PRIVATE galois doctest::doctest)
target_compile_options(
  rsiscool-tests
  PRIVATE -Og
          -g
          -Wall
          -Wextra
          -Wshadow
          -pedantic
          -Wno-missing-field-initializers
          -Wno-deprecated-declarations)

