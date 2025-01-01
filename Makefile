.PHONY: build

EXECUTABLE = test/rsiscool-tests

tests:
	mkdir -p test/
	cd test && cmake .. -DCMAKE_TOOLCHAIN_FILE=clang.cmake
	cmake --build test --target rsiscool-tests

	cp test/compile_commands.json .

module:
	sh build_module.sh

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
