.PHONY: build

EXECUTABLE = test/rsiscool-tests

all:
	mkdir -p build/
	cd build && emcmake cmake ..
	cmake --build build --target rsiscool

	mkdir -p test/
	cd test && cmake .. -DCMAKE_TOOLCHAIN_FILE=clang.cmake
	cmake --build test --target rsiscool-tests

	mkdir -p build-workers/
	cd build-workers && emcmake cmake ..
	cmake --build build-workers --target rsiscool_workers

	cp test/compile_commands.json .

	cp build/rsiscool.js dist/
	cp build/rsiscool.d.ts dist/
	cp build/rsiscool.wasm dist/

	cp build-workers/rsiscool_workers.js dist/
	cp build-workers/rsiscool_workers.d.ts dist/

build:
	mkdir -p test/
	cd test && cmake .. -DCMAKE_TOOLCHAIN_FILE=clang.cmake
	cmake --build test --target rsiscool-tests

	cp test/compile_commands.json .

run:
	$(EXECUTABLE)
