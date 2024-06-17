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

	# do this properly
	cp build/rsiscool.js ../jsQR/src/wasm/rsiscool.js
	cp build/rsiscool.wasm ../jsQR/src/wasm/rsiscool.wasm

build:
	mkdir -p test/
	cd test && cmake .. -DCMAKE_TOOLCHAIN_FILE=clang.cmake
	cmake --build test --target rsiscool-tests

	cp test/compile_commands.json .

run:
	$(EXECUTABLE)
