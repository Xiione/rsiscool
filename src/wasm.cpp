#include <vector>
#include <emscripten/bind.h>
#include "ReedSolomon.hpp"

struct DecodeResult {
  int errors;
  std::optional<std::vector<uint8_t>> bytesCorrected;
};

// bytes: number[]
DecodeResult decodeWASM(const emscripten::val &bytes, int twoS) {
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
    return {res, std::nullopt};
  }

  return {res, v};
}

EMSCRIPTEN_BINDINGS(rsiscool) {
  emscripten::value_object<DecodeResult>("DecodeResult")
    .field("errors", &DecodeResult::errors)
    .field("bytesCorrected", &DecodeResult::bytesCorrected);
  emscripten::register_vector<uint8_t>("Uint8Vector");
  emscripten::register_optional<std::vector<uint8_t>>();
  emscripten::function("decodeWASM", &decodeWASM);
}
