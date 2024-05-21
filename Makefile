.PHONY: build

# TYPE = release
TYPE = debug

build:
	cmake --build build/$(TYPE) --target decoder
	cp build/$(TYPE)/compile_commands.json .

run:
	./build/$(TYPE)/qriscool
