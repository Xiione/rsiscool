#pragma once

#include <optional>
#include <vector>

#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <galois/GaloisField.h>

constexpr uint ALPHA = 2;

struct ReedSolomon {
  uint Q;
  uint N;
  uint K;
  // galois field with LUT used for mul and exp
  galois::GaloisField gfLUT;
  uint primPowOffset;

  NTL::GF2EX genPoly;

  ReedSolomon(uint Q, uint K, const std::vector<uint> &primPoly,
              uint primPowOffset = 0);

  NTL::GF2EX encode(std::vector<NTL::GF2E> &u);

  // returns nullopt if at least N - K errors => error cannot be solved for
  std::optional<NTL::GF2EX> decodePGZ(const NTL::GF2EX &r,
                                      size_t const *resErrs = nullptr);

  std::optional<NTL::GF2EX> decodeBM(const NTL::GF2EX &r,
                                     size_t const *resErrs = nullptr);
};

NTL::GF2X intToGF2X(uint x);
NTL::GF2E intToGF2E(uint x);
NTL::GF2X vecToGF2X(const std::vector<uint> v);
