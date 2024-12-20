#pragma once

#include <optional>
#include <vector>

#include <galois/GaloisField.h>
#include <galois/GaloisFieldElement.h>
#include <galois/GaloisFieldPolynomial.h>

constexpr int ALPHA = 2;
const std::vector<int> SYMBOL_PRIMPOLY = {1, 0, 1, 1, 1, 0, 0, 0, 1};

struct ReedSolomon {
  int N;
  int K;
  int primRootOffset;

#ifdef RSISCOOL_ENCODE
  galois::GaloisFieldPolynomial genPoly;
#endif

  ReedSolomon(int N, int K, int primPowOffset = 0);

  // #ifdef RSISCOOL_ENCODE
  galois::GaloisFieldPolynomial
  encode(std::vector<galois::GaloisFieldElement> &input);
  // #endif

  // returns nullopt if at least N - K errors => error cannot be solved for
  // uses linear solving for error vals
  // std::optional<galois::GaloisFieldPolynomial>
  // decodePGZ(galois::GaloisFieldElement received, int *const resErrs =
  // nullptr);

  // uses forney algorithm
  std::optional<galois::GaloisFieldPolynomial>
  decodeBM(galois::GaloisFieldPolynomial received,
           int *const resErrs = nullptr);

  // chien search from wikipedia
  // give coefficients with index 0 as coefficient of x^0 term
  std::vector<int> findRootPows(std::vector<galois::GaloisFieldElement> coeffs);
  // give polynomial directly
  std::vector<int> findRootPows(const galois::GaloisFieldPolynomial &loc);

  // std::optional<std::vector<galois::GaloisFieldElement>>
  // solveErrorVals(const std::vector<galois::GaloisFieldElement> &syndromes,
  //                const std::vector<int> &errorLocs);

  std::optional<std::vector<galois::GaloisFieldElement>>
  solveErrorValsForney(galois::GaloisFieldPolynomial &locator,
                       std::vector<galois::GaloisFieldElement> &syndromes,
                       const std::vector<int> &errorLocs);

  // returns whether not poly evaluated at roots of generator polynomial are all
  // zero
  bool checkSyndromes(galois::GaloisFieldPolynomial poly);
};

galois::GaloisFieldPolynomial bytesToGF2EX(std::vector<uint8_t> &bytes);
int decodeBytes(std::vector<uint8_t> &bytes, int twoS);
bool validateBytes(std::vector<uint8_t> &bytes, int twoS);
