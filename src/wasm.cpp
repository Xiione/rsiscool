#include "ReedSolomon.hpp"
#include "wasm.hpp"
#include <emscripten/bind.h>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

DecodeResult decodeWASM(const std::string &bytes, int twoS) {
  ReedSolomon rs(bytes.size(), bytes.size() - twoS);
  std::vector<uint8_t> v(bytes.begin(), bytes.end());

  auto [errors, locs, vals] = rs.decodeBM(bytesToGF2EX(v));

  if (!locs || !vals) {
    return {errors, std::nullopt};
  }

  for (int i = 0; i < locs->size(); ++i) {
    // reverse back to jsqr order
    v[rs.N - locs->at(i) - 1] ^= vals->at(i).poly();
  }

  return {errors, v};
}

bool validateWASM(const std::string &bytes, int twoS) {
  ReedSolomon rs(bytes.size(), bytes.size() - twoS);
  return rs.checkSyndromes(bytesToGF2EX(std::span<const uint8_t>(
      reinterpret_cast<const uint8_t *>(bytes.data()), bytes.size())));
}

#ifndef RSISCOOL_NOEMBIND
EMSCRIPTEN_BINDINGS(rsiscool) {
  emscripten::value_object<DecodeResult>("DecodeResult")
      .field("errors", &DecodeResult::errors)
      .field("bytesCorrected", &DecodeResult::bytesCorrected);
  emscripten::register_vector<uint8_t>("Uint8Vector");
  emscripten::register_optional<std::vector<uint8_t>>();
  emscripten::function("decodeWASM", &decodeWASM);
  emscripten::function("validateWASM", &validateWASM);
}
#endif
