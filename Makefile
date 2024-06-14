.PHONY: build

EXECUTABLE = build/rsiscool-test

build:
	cd build && emcmake cmake ..
	cmake --build build --target rsiscool

	cd test && cmake ..
	cmake --build test --target rsiscool-tests

	cp test/compile_commands.json .

run:
	$(EXECUTABLE)
