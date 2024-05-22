cmake_minimum_required(VERSION 3.5...3.28)

project(
  qriscool
  VERSION 0.0.1
  DESCRIPTION "QR Is Cool by Hamilton Wang")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # clangd

add_subdirectory(deps/galois)
add_subdirectory(deps/eigen)

set(APP_SRC 
    src/main.cpp
    src/ReedSolomon.cpp
)

add_executable(qriscool ${APP_SRC})

target_link_libraries(qriscool PRIVATE galois eigen)
