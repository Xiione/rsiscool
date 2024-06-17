For posterity:
Versions using NTL/GF2X were installed by cloning the repo sources into deps/ and then running emconfigure ... HOST_CC=gcc. NTL fails to produce mach_desc.h so it's copied from the given (??) mach_desc.win file since it mostly matches the output of configure.
The brew versions of the libraries were also installed because NTL performs some platform-specific config at installation time which differs between when you want to compile the doctest executable and the wasm library.
