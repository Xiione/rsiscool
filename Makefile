.PHONY: build

# TYPE = release
TYPE = debug

install:
	cd deps/givaro

	# also add --with-gmp=<path-to-gmp>/x.x.x/ if gmp not found
	./autogen.sh --prefix=$(shell pwd) --disable-shared ABI=64 CFLAGS="-m64 -O2"

build:
	cmake --build build/$(TYPE) --target qriscool
	cp build/$(TYPE)/compile_commands.json .

run:
	./build/$(TYPE)/qriscool
