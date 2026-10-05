.PHONY: build

EXECUTABLE = test/rsiscool-tests
EMSCRIPTEN_ROOT := $(shell em-config EMSCRIPTEN_ROOT)
EMSCRIPTEN_BUILD_ID := $(notdir $(patsubst %/,%,$(dir $(EMSCRIPTEN_ROOT))))
BUILD_DIR := build/$(EMSCRIPTEN_BUILD_ID)
WORKERS_BUILD_DIR := build-workers/$(EMSCRIPTEN_BUILD_ID)

tests:
	mkdir -p test/
	cmake -B test -S . -DCMAKE_TOOLCHAIN_FILE=clang.cmake
	cmake --build test --target rsiscool-tests

	cp test/compile_commands.json .

module:
	mkdir -p $(BUILD_DIR)
	emcmake cmake -B $(BUILD_DIR) -S .
	cmake --build $(BUILD_DIR) --target rsiscool

	cp $(BUILD_DIR)/rsiscool.js dist/
	cp $(BUILD_DIR)/rsiscool.d.ts dist/
	cp $(BUILD_DIR)/rsiscool.wasm dist/

	mkdir -p $(WORKERS_BUILD_DIR)
	emcmake cmake -B $(WORKERS_BUILD_DIR) -S .
	cmake --build $(WORKERS_BUILD_DIR) --target rsiscool_workers

	cp $(WORKERS_BUILD_DIR)/rsiscool_workers.js dist/
	cp $(WORKERS_BUILD_DIR)/rsiscool_workers.d.ts dist/
	cp $(WORKERS_BUILD_DIR)/rsiscool_workers.wasm dist/

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
