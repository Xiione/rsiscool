#pragma once

#include <optional>
#include <vector>

#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/matrix.h>

constexpr int ALPHA = 2;
const std::vector<int> SYMBOL_PRIMPOLY = {1, 0, 1, 1, 1, 0, 0, 0, 1};

struct ReedSolomon {
  int N;
  int K;
  int primRootOffset;

  NTL::GF2EX genPoly;

  ReedSolomon(int N, int K, int primPowOffset = 0);

#ifdef RSISCOOL_ENCODE
  NTL::GF2EX encode(std::vector<NTL::GF2E> &input);
#endif

  // returns nullopt if at least N - K errors => error cannot be solved for
  // uses linear solving for error vals
  std::optional<NTL::GF2EX> decodePGZ(NTL::GF2EX received,
                                      int *const resErrs = nullptr);

  // uses forney algorithm
  std::optional<NTL::GF2EX> decodeBM(NTL::GF2EX received,
                                     int *const resErrs = nullptr);

  // chien search from wikipedia
  // give coefficients with index 0 as coefficient of x^0 term
  std::vector<int> findRootPows(const NTL::Vec<NTL::GF2E> &coeffs);

  std::optional<std::vector<NTL::GF2E>>
  solveErrorVals(const std::vector<NTL::GF2E> &syndromes,
                 const std::vector<int> &errorLocs);

  std::optional<std::vector<NTL::GF2E>>
  solveErrorValsForney(const NTL::GF2EX &locator,
                       const std::vector<NTL::GF2E> &syndromes,
                       const std::vector<int> &errorLocs);
};

void initGF2E();
std::optional<std::vector<uint8_t>>
decodeBytes(const std::vector<uint8_t> &bytes, int twoS);
