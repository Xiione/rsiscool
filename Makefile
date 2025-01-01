.PHONY: build

EXECUTABLE = test/rsiscool-tests

tests:
	mkdir -p test/
	cd test && cmake .. -DCMAKE_TOOLCHAIN_FILE=clang.cmake
	cmake --build test --target rsiscool-tests

	cp test/compile_commands.json .

module:
	mkdir -p build/
	emcmake cmake -B build -S .
	cmake --build build --target rsiscool

	cp build/rsiscool.js dist/
	cp build/rsiscool.d.ts dist/
	cp build/rsiscool.wasm dist/

	mkdir -p build-workers/
	emcmake cmake -B build-workers -S .
	cmake --build build-workers --target rsiscool_workers

	cp build-workers/rsiscool_workers.js dist/
	cp build-workers/rsiscool_workers.d.ts dist/
	ln -f build/rsiscool.wasm dist/rsiscool_workers.wasm

build:
	make tests

all:
	make tests
	make module

clean:
	rm -rf ./build
	rm -rf ./build-workers
	rm -rf ./test

run:
	$(EXECUTABLE)
