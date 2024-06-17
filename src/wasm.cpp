#include <algorithm>
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


  auto res = decodeBytes(v, twoS);
  if (!res) {
    return emscripten::val::null();
  }

  emscripten::val arr =
      emscripten::val::global("Uint8ClampedArray").new_(res->size());
  arr.call<void>("set", emscripten::val(emscripten::typed_memory_view(
                            res->size(), res->data())));

  return arr;
}

EMSCRIPTEN_BINDINGS(rsiscool) {
  emscripten::function("initGF2E", &initGF2E);
  emscripten::function("decodeWASM", &decodeWASM);
}
