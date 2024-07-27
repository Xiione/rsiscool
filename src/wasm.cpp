#include "ReedSolomon.hpp"
#include <emscripten/bind.h>
#include <vector>

struct DecodeResult {
  int errors;
  std::optional<std::vector<uint8_t>> bytesCorrected;
};

// bytes: number[]
DecodeResult decodeWASM(const emscripten::val &bytes, int twoS) {
  // https://github.com/emscripten-core/emscripten/issues/5519#issuecomment-333302296
  emscripten::val memory = emscripten::val::module_property("HEAPU8")["buffer"];

  std::vector<uint8_t> v =
      emscripten::convertJSArrayToNumberVector<uint8_t>(bytes);

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
