.PHONY: build

EXECUTABLE = test/rsiscool-tests

all:
	mkdir -p build/
	cd build && emcmake cmake ..
	cmake --build build --target rsiscool

	mkdir -p test/
	cd test && cmake .. -DCMAKE_TOOLCHAIN_FILE=clang.cmake
	cmake --build test --target rsiscool-tests

	cp test/compile_commands.json .

	cp build/rsiscool.js .
	cp build/rsiscool.wasm .

build:
	mkdir -p test/
	cd test && cmake .. -DCMAKE_TOOLCHAIN_FILE=clang.cmake
	cmake --build test --target rsiscool-tests

	cp test/compile_commands.json .

run:
	$(EXECUTABLE)
