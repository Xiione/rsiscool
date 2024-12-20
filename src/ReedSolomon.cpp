#include <algorithm>
#include <cassert>
#include <numeric>
#include <optional>
#include <sys/types.h>
#include <vector>

#include <galois/GaloisField.h>
#include <galois/GaloisFieldElement.h>
#include <galois/GaloisFieldPolynomial.h>
#include <galois/GlobalGaloisField.h>

#include "ReedSolomon.hpp"

using GF2 = galois::GaloisField;
using GF2E = galois::GaloisFieldElement;
using GF2EX = galois::GaloisFieldPolynomial;

const uint primPoly[] = {1, 0, 1, 1, 1, 0, 0, 0, 1};
GF2 GF(8, primPoly);

inline GF2E coeff(const GF2EX &a, long i) {
  if (i > a.deg())
    return GF2E(0);
  return a[i];
}

GF2EX bytesToGF2EX(std::vector<uint8_t> &bytes) {
  std::vector<GF2E> coeffs(bytes.size());
  // do reversing here
  std::transform(bytes.rbegin(), bytes.rend(), coeffs.begin(),
                 [](uint8_t b) { return GF2E(b); });

  return GF2EX(bytes.size() - 1, coeffs.data());
}

int decodeBytes(std::vector<uint8_t> &bytes, int twoS) {
  ReedSolomon rs(bytes.size(), bytes.size() - twoS);

  int errors = -1;
  auto res = rs.decodeBM(bytesToGF2EX(bytes), &errors);

  if (res) {
    for (int i = 0; i < rs.N; ++i) {
      bytes[rs.N - i - 1] = coeff(res.value(), i).poly();
    }
  }
  return errors;
}

bool validateBytes(std::vector<uint8_t> &bytes, int twoS) {
  ReedSolomon rs(bytes.size(), bytes.size() - twoS);

  GF2EX code = bytesToGF2EX(bytes);
  return rs.checkSyndromes(code);
}

ReedSolomon::ReedSolomon(int N, int K, int offset)
    : N(N), K(K), primRootOffset(offset) {
  // init generator polynomial
#ifdef RSISCOOL_ENCODE
  int genDeg = N - K;
  std::vector<GF2E> coeffInit = {GF2E(1)};
  genPoly = GF2EX(0, coeffInit.data());
  // (z - a^i)
  std::vector<GF2E> minPolyCoeff = {GF2E(0), GF2E(1)};

  for (int i = 0; i < genDeg; ++i) {
    minPolyCoeff[0] = GF->alpha(offset + i);
    // *= (z - a^(i + 1))
    genPoly *= GF2EX(&gf, 1, minPolyCoeff.data());
  }
#endif
}

#ifdef RSISCOOL_ENCODE
GF2EX ReedSolomon::encode(std::vector<GF2E> &input) {
  GF2EX c(input.size(), input.data());

  // u(z) * z^(n - k)
  c <<= N - K;
  // u(z) * z^(n - k) - R_{g(z)}[u(z) * z^(n-k)]
  c -= c % genPoly;
  return c;
}
#endif

std::optional<GF2EX> ReedSolomon::decodeBM(GF2EX rec, int *const resErrs) {
  // maximum correctable errors
  unsigned t = (N - K + 1) / 2;
  std::vector<GF2E> syndromes(N - K);

  for (int i = 0; i < N - K; ++i) {
    syndromes[i] = rec(GF.alpha(primRootOffset + i));
  }

  GF2E one(1);
  GF2EX loc_m(0, &one), loc_r(0, &one);
  int m = -1;
  int len_m = 0, len_r = 0;
  GF2E delta_m(1);
  for (int r = 0; r < (int)(N - K); ++r) {
    GF2E delta_r(0);
    for (int h = 0; h <= len_r; ++h) {
      delta_r += coeff(loc_r, h) * syndromes[r - h];
    }

    if (delta_r == 0) {
      // no change, current solution still still solves for (r + 1)th syndrome
      continue;
    }

    GF2EX locTmp = loc_r - ((delta_r * delta_m.inverse()) * (loc_m << (r - m)));
    int lenTmp = std::max(len_r, len_m + r - m);

    // new maximum m - l_m prev solution
    if (r - len_r > m - len_m) {
      loc_m = loc_r;
      m = r;
      len_m = len_r;
      delta_m = delta_r;
    }
    len_r = lenTmp;
    loc_r = locTmp;
  }

  // too many errors to correct
  if (loc_r.deg() > t) {
    if (resErrs != nullptr)
      *resErrs = -1;
    return std::nullopt;
  }

  std::vector<int> locs = findRootPows(loc_r);

  auto vals = solveErrorValsForney(loc_r, syndromes, locs);

  if (!vals) {
    if (resErrs != nullptr)
      *resErrs = -1;
    return std::nullopt;
  }

  for (int i = 0; i < locs.size(); ++i) {
    if (locs[i] > rec.deg())
      rec.set_degree(locs[i]);
    rec[locs[i]] -= vals.value()[i];
  }

  // final sanity check, are syndromes zero now?
  if (!checkSyndromes(rec)) {
    if (resErrs != nullptr)
      *resErrs = -1;
    return std::nullopt;
  }

  if (resErrs != nullptr)
    *resErrs = locs.size();
  return rec;
}

std::vector<int> ReedSolomon::findRootPows(std::vector<GF2E> v) {
  std::vector<int> res;
  for (int i = 0; i < N; ++i) {
    if (std::accumulate(v.begin(), v.end(), GF2E(0)) == 0) {
      res.push_back(i);
    }

    for (int j = 0; j < v.size(); ++j) {
      v[j] *= GF.alpha(j);
    }
  }
  return res;
}

std::vector<int> ReedSolomon::findRootPows(const GF2EX &loc) {
  std::vector<GF2E> v(loc.deg() + 1);
  for (int i = 0; i <= loc.deg(); ++i)
    v[i] = loc[loc.deg() - i];
  std::vector<int> res;
  for (int i = 0; i < N; ++i) {
    if (std::accumulate(v.begin(), v.end(), GF2E(0)) == 0) {
      res.push_back(i);
    }

    for (int j = 0; j < v.size(); ++j) {
      v[j] *= GF.alpha(j);
    }
  }
  return res;
}

std::optional<std::vector<GF2E>>
ReedSolomon::solveErrorValsForney(GF2EX &locator, std::vector<GF2E> &syndromes,
                                  const std::vector<int> &errorLocs) {
  // initialize with syndrome polynomial
  GF2EX O(syndromes.size() - 1, syndromes.data());

  O = O * locator;
  O.set_degree(N - K - 1);
  O.simplify();

  GF2EX locatorDiff = locator.derivative();

  std::vector<GF2E> vals;
  vals.resize(errorLocs.size());
  for (int i = 0; i < errorLocs.size(); ++i) {
    GF2E loc = GF2E(GF.alpha(errorLocs[i]));
    GF2E locInv = GF2E(loc.inverse());

    GF2E l = locatorDiff(locInv);
    if (l == 0)
      return std::nullopt;

    // Gill, John EE387, source from wikipedia. for non-1 offsets the first
    // term is necessary
    if (primRootOffset == 0)
      vals[i] = loc * O(locInv) / l;
    else
      vals[i] = GF.exp(loc.poly(), 1 - (int)primRootOffset) * O(locInv) / l;
  }

  return vals;
}

bool ReedSolomon::checkSyndromes(GF2EX poly) {
  for (int i = 0; i < N - K; ++i) {
    if (poly(GF.alpha(primRootOffset + i)) != 0) {
      return false;
    }
  }
  return true;
}
