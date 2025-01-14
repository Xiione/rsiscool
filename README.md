# rsiscool
rsiscool is a small library implementing several algorithms related to [Reed-Solomon codes](https://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction), including the Berlekamp-Massey algorithm for decoding Reed-Solomon codes. It is compiled to WebAssembly using [Emscripten](https://emscripten.org/) and uses [Embind](https://emscripten.org/docs/porting/connecting_cpp_and_javascript/embind.html) to create function bindings for calling from Javascript code. The library is configured to precompute lookup tables for only the finite field GF(2⁸), as the QR code standard exclusively utilizes Reed-Solomon codes of this field. 

To see rsiscool in action, visit [qris.cool](https://qris.cool)! This library is used in the client-side application as well as in the server-side validation code.

## Function bindings
Function bindings are registered in the `EMSCRIPTEN_BINDINGS` macro invocation in `src/wasm.cpp`.

### `DecodeResult decodeWASM(const std::string &bytes, int twoS)`
`bytes`: An `ArrayBuffer`, `Uint8Array`, `Uint8ClampedArray`, `Int8Array` of length **N** representing the GF(2⁸) coefficients of an (**N**, K) Reed-Solomon code, in _reverse_ order of degree, i.e. `bytes[0]` is the coefficient of x^(**N** - 1) while `bytes[N - 1]` is the coefficient of x^0.

`twoS`: An integer `number` representing the number of redundancy symbols in the code, in other words the value N - K for the (N, K) Reed-Solomon code. The naming is related to the fact that the minimum Hamming distance of a code is equal to about half the degree of the generating polynomial of the code, and thus the count of redundancy symbols.

Returns a `struct DecodeResult` of two elements. One, an `optional<vector<uint8_t>>` `bytesCorrected` containing a vector of the coefficients of the corrected Reed-Solomon code (in the same order as the input) if correction was successful, otherwise `std::nullopt`. Two, an integer `errors`, the count of errors corrected. This value is undefined if correction of the input code fails.

The `vector<uint8_t>` class is registered via Embind and the member functions `get`, `set`, `size`, `push_back`, `resize`, as well as an additional cleanup function `delete` calling the destructor internally are exposed and callable from Javascript on the Javascript object.

### `bool validateWASM(const std::string &bytes, int twoS)`
`bytes`: Same as above. See `decodeWASM`.  
`twoS`: Same as above. See `decodeWASM`.

Returns true if the given code is a valid Reed-Solomon code. The code is valid if and only if the generator polynomial's roots are also roots of the code polynomial, thus the given code is evaluated at these points and only if all the outputs are 0 does this function return true.

## CMake build targets
The tests and two WASM outputs are each separated to prevent CMake from mixing platform-specific object output. See [CMakeLists.txt](https://github.com/Xiione/rsiscool/blob/main/CMakeLists.txt).
The WASM outputs are tracked in this repo because Cloudflare's build environment for Pages and Workers seems to have the Emscripten toolchain installed improperly, and also because the build times are slow.


### rsiscool-tests (Native)
The executable built by this target runs a small set of regression tests on the two Embinded functions. The macro `RSISCOOL_NOEMBIND` is defined to exclude the `EMSCRIPTEN_BINDINGS` invocation in this build. Uses [doctest](https://github.com/doctest/doctest).

### rsiscool (WASM)
This build is used in the Svelte client. In the interest of size the WASM binary is included separately from the JS glue code emitted by em++, Vite is able to properly bundle the binary to the expected location.  
Link flags:
- `--no-entry`: This build only emits function bindings.
- `SINGLE_FILE=0`: The WASM should not be embedded as a base64 string.
- `EXPORT_ES6=1`, `MODULARIZE=1`: Modules are nice to work with.
- `EMBIND_STD_STRING_IS_UTF8=0`: Do not assume Embinded functions with std::string parameters are in UTF-8 encoding as we do pass arbitrary `Uint8Array`s.
- `--emit-tsd 'rsiscool.d.ts'`: Makes development slightly easier with types.
- `-lembind`: We need Embind linked.

### rsiscool_workers (WASM)
This build is used in the Cloudflare Workers code used to [validate user-submitted QR codes](https://qris.cool/privacy).  
Link flags:
- `--no-entry`: Same as `rsiscool`.
- `SINGLE_FILE=0`: The WASM should not be embedded as a base64 string. For security reasons, [Cloudflare actually disallows it explicitly](https://github.com/thx/resvg-js/issues/307) even though it works locally as expected using Wrangler. The solution is to import the WASM binary as a module in the worker code ([see this example](https://developers.cloudflare.com/workers/runtime-apis/webassembly/javascript/#use-from-javascript)) and [instantiate the WASM module manually](https://github.com/Xiione/jsQR/blob/master/src/decoder/reedsolomon/index.ts#L14). We MUST use the factory function exported from the glue code because the `instantiateWasm` callback passes the WASM instantiation several crucial imports that the glue code normally handles when it locates and instantiates the WASM on its own. 
- `DYNAMIC_EXECUTION=0`: Disables the emitting of glue JS utilizing dynamic execution of arbitrary code. Similar to `SINGLE_FILE=1`, Cloudflare also disallows this but here Wrangler does actually throw an error when you use it in local development. 
- `EMBIND_AOT=1`: According to Emscripten docs, this provides a speed boost to Embind JS invoker functions by generating them at compile time, combined with `DYNAMIC_EXECUTION=0`, we achieve performance identical to the client WASM module.
- `ENVIRONMENT='web'`: This prevents the module from emitting glue JS meant for compatibility with multiple environments (which would ultimately cause initialization to error, like `import.meta.url` being undefined) Even though Workers is a node-like environment, `'web'` is what works.
- `EXPORT_ES6=1`, `MODULARIZE=1`: Same as `rsiscool`.
- `EMBIND_STD_STRING_IS_UTF8=0`: Same as `rsiscool`.
- `--emit-tsd 'rsiscool_workers.d.ts'`: The emitted types are identical to `rsiscool.d.ts` but we take a union of the two modules' types to make the language server/compiler happy.
- `-lembind`: Same as `rsiscool`.

## Credits
- The C++ library implementating Galois fields and Galois field polynomials included directly in this repository was made by Arash Partow. It was downloaded from the ["Galois Field Arithmetic Library" webpage](https://www.partow.net/projects/galois/) on Arash Partow's website.  
- [doctest](https://github.com/doctest/doctest) is used for the tests in this project. For local development it was installed via homebrew.  
- Credit to a minimal example given in [a repo created by robertaboukhalil](https://github.com/robertaboukhalil/cf-workers-emscripten) for helping me get WASM to work in Workers.
- Prior versions of this library also used [NTL](https://github.com/libntl/ntl) or [Givaro](https://github.com/linbox-team/givaro) for their implementations of finite fields and linear algebra in fields.

A big thank you to each of these developers/groups for their work!
