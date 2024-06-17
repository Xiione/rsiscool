#include <vector>
#include <emscripten/bind.h>
#include "ReedSolomon.hpp"

// bytes: number[]
emscripten::val decodeWASM(const emscripten::val &bytes, int twoS) {
  if (bytes == emscripten::val::null()) {
    return emscripten::val::null();
  }

  // https://github.com/emscripten-core/emscripten/issues/5519#issuecomment-333302296
  emscripten::val memory = emscripten::val::module_property("HEAPU8")["buffer"];
  unsigned int length = bytes["length"].as<unsigned int>();

  std::vector<uint8_t> v;
  v.resize(length);

  emscripten::val inMemView =
      emscripten::val::global("Uint8Array")
          .new_(memory, reinterpret_cast<uintptr_t>(v.data()), length);
  inMemView.call<void>("set", bytes);


  int res = decodeBytes(v, twoS);
  if (res == -1) {
    return emscripten::val::null();
  }

  emscripten::val arr =
      emscripten::val::global("Uint8ClampedArray").new_(v.size());
  arr.call<void>("set", emscripten::val(emscripten::typed_memory_view(
                            v.size(), v.data())));

  return arr;
}

EMSCRIPTEN_BINDINGS(rsiscool) {
  emscripten::function("decodeWASM", &decodeWASM);
}
