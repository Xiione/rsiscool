#pragma once

#include <optional>
#include <string>
#include <vector>

struct DecodeResult {
  size_t errors;
  std::optional<std::vector<uint8_t>> bytesCorrected;
};

DecodeResult decodeWASM(const std::string &bytes, int twoS);
bool validateWASM(const std::string &bytes, int twoS);
