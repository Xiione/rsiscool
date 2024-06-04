.PHONY: build

# TYPE = release
TYPE = debug
EXECUTABLE = ./build/debug/qriscool


# install:
 	# pre-install:
	# brew install parkerdiamond/gf2x/gf2x
	# brew install ntl
	# brew edit ntl and add option NTL_GF2X_LIB=on
	# brew reinstall --build-from-source ntl 

build:
	cmake --build build/$(TYPE) --target qriscool
	cp build/$(TYPE)/compile_commands.json .

run:
	$(EXECUTABLE)
