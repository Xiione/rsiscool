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
  uint primPowOffset;
  std::vector<NTL::GF2E> primPow;

  NTL::GF2EX genPoly;

  ReedSolomon(uint Q, uint K, uint primPowOffset = 0);

  NTL::GF2EX encode(std::vector<NTL::GF2E> &u);

  // returns nullopt if at least N - K errors => error cannot be solved for
  std::optional<NTL::GF2EX> decodePGZ(NTL::GF2EX r,
                                      uint *const resErrs = nullptr);

  std::optional<NTL::GF2EX> decodeBM(NTL::GF2EX r,
                                     uint *const resErrs = nullptr);

  // chien search from wikipedia
  // give coefficients with index 0 as coefficient of x^0 term
  std::vector<uint> findRoots(NTL::Vec<NTL::GF2E> v);
};

NTL::GF2X intToGF2X(uint x);
NTL::GF2E intToGF2E(uint x);
NTL::GF2X vecToGF2X(const std::vector<uint> v);
uint GF2EtoInt(const NTL::GF2E &x);
