#pragma once

#include <optional>
#include <vector>

#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/matrix.h>
#include <galois/GaloisField.h>

constexpr uint ALPHA = 2;

struct ReedSolomon {
  uint Q;
  uint N;
  uint K;
  uint primRootOffset;
  std::vector<NTL::GF2E> primPow;

  NTL::GF2EX genPoly;

  ReedSolomon(uint Q, uint K, uint primPowOffset = 0);

  NTL::GF2EX encode(std::vector<NTL::GF2E> &input);

  // returns nullopt if at least N - K errors => error cannot be solved for
  // uses linear solving for error vals
  std::optional<NTL::GF2EX> decodePGZ(NTL::GF2EX received,
                                      int *const resErrs = nullptr);

  // uses forney algorithm
  std::optional<NTL::GF2EX> decodeBM(NTL::GF2EX received,
                                     int *const resErrs = nullptr);

  // chien search from wikipedia
  // give coefficients with index 0 as coefficient of x^0 term
  std::vector<uint> findRootPows(NTL::Vec<NTL::GF2E> coeffs);

  std::optional<std::vector<NTL::GF2E>>
  solveErrorVals(const std::vector<NTL::GF2E> &syndromes,
                 const std::vector<uint> &errorLocs);

  std::optional<std::vector<NTL::GF2E>>
  solveErrorValsForney(const NTL::GF2EX &locator,
                       const std::vector<NTL::GF2E> &syndromes,
                       const std::vector<uint> &errorLocs);
};

NTL::GF2X intToGF2X(uint x);
NTL::GF2E intToGF2E(uint x);
NTL::GF2X vecToGF2X(const std::vector<uint> v);
uint GF2EtoInt(const NTL::GF2E &x);
