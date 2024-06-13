# install:


build:
	emcmake cmake . -B dist
	emmake make -C dist

	# cmake --build build/$(TYPE) --target rsiscool
	cp dist/compile_commands.json .
