#include <algorithm>
#include <cassert>
#include <numeric>
#include <optional>
#include <unordered_set>
#include <vector>

#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2E.h>

#include "ReedSolomon.hpp"

// #define getRep(x) x._GF2E__rep.xrep.rep
// #define GF2EtoInt(x) (getRep(x) ? *getRep(x) : 0)

std::vector<NTL::GF2E> primPow;

// random
inline NTL::GF2X bytetoGF2X(uint8_t x) { return NTL::GF2XFromBytes(&x, 1); }

inline NTL::GF2E byteToGF2E(uint8_t x) { return NTL::to_GF2E(bytetoGF2X(x)); }

NTL::GF2X vecToGF2X(const std::vector<int> v) {
  size_t numBytes = (v.size() + 7) / 8;
  std::vector<uint8_t> res(numBytes, 0);

  for (size_t i = 0; i < v.size(); ++i) {
    res[i / 8] |= (v[i] != 0) << (i % 8);
  }
  return NTL::GF2XFromBytes(res.data(), numBytes);
}

// evil random error where gf2e value changes when passed as an argument
// we will use a macro instead
// inline uint8_t GF2EtoInt(const NTL::GF2E &x) {
//   NTL::GF2X f = x._GF2E__rep;
//   if (f.xrep.rep)
//     return *f.xrep.rep;
//   return 0;
// }

void reduce(NTL::Mat<NTL::GF2E> &A, int rows) {
  for (int i = rows - 1; i >= 0; --i) {
    A[i] *= NTL::inv(A[i][i]);

    for (int j = 0; j < i; ++j) {
      A[j] -= A[j][i] * A[i];
    }
  }
}

// check if augmented matrix representing a system is independent and consistent
bool checkSystem(NTL::Mat<NTL::GF2E> &A, int rank) {
  // check that each equation is consistent
  for (int i = 0; i < A.NumRows(); ++i) {
    bool allZeroCoeffs = true;
    for (int j = 0; j < rank; ++j) {
      allZeroCoeffs = allZeroCoeffs && NTL::IsZero(A[i][j]);
    }
    // system was not reduced properly, NTL::gauss or my reduce's fault?
    if (i >= rank && !(allZeroCoeffs && NTL::IsZero(A[i][rank])))
      return false;
    // all zero coeffs but right side is nonzero -> inconsistent
    if (allZeroCoeffs && !NTL::IsZero(A[i][rank]))
      return false;
  }
  return true;
}

void initGF2E() {
  static bool initialized = false;
  if (initialized)
    return;
  initialized = true;

  NTL::GF2X primPoly = vecToGF2X(SYMBOL_PRIMPOLY);
  std::cerr << "[rsiscool] Using primitive polynomial: " << primPoly
            << std::endl;
  NTL::GF2E::init(primPoly);

  NTL::GF2E a = byteToGF2E(ALPHA);
  NTL::GF2E x = NTL::GF2E(1);

  primPow.resize(256);
  for (int i = 0; i < 256; ++i) {
    primPow[i] = x;
    x *= a;
  }
}

int decodeBytes(std::vector<uint8_t> &bytes, int twoS) {
  ReedSolomon rs(bytes.size(), bytes.size() - twoS);
  NTL::GF2EX rec;
  for (int i = 0; i < rs.N; ++i) {
    // do reversing here
    NTL::SetCoeff(rec, i, byteToGF2E(bytes[rs.N - i - 1]));
  }

  int errors = -2;
  auto res = rs.decodeBM(rec, &errors);
  if (!res)
    return errors;

  for (int i = 0; i < rs.N; ++i) {
    auto f = NTL::coeff(*res, i)._GF2E__rep.xrep;
    bytes[rs.N - i - 1] = f.rep ? (uint8_t)*f.rep : 0;
  }
  return errors;
}

ReedSolomon::ReedSolomon(int N, int K, int offset)
    : N(N), K(K), primRootOffset(offset) {
  // deduce extension field power
  // assert(Q == 1 << (int)log2(Q));

  // init generator polynomial
#ifdef RSISCOOL_ENCODE
  int genDeg = N - K;
  NTL::GF2EX primMinPoly;
  // (z - a^i)
  NTL::Vec<NTL::GF2E> roots;
  roots.SetLength(genDeg);
  for (size_t i = 0; i < genDeg; ++i) {
    roots[i] = primPow[offset + i];
  }

  genPoly = NTL::BuildFromRoots(roots);
  assert(NTL::deg(genPoly) == genDeg);
#endif
}

#ifdef RSISCOOL_ENCODE
NTL::GF2EX ReedSolomon::encode(std::vector<NTL::GF2E> &input) {
  NTL::GF2EX c;
  for (int i = 0; i < input.size(); ++i)
    NTL::SetCoeff(c, i, input[i]);

  // u(z) * z^(n - k)
  c <<= N - K;
  assert(NTL::deg(c) == N - 1);

  // u(z) * z^(n - k) - R_{g(z)}[u(z) * z^(n-k)]
  c -= c % genPoly;
  return c;
}
#endif

std::optional<NTL::GF2EX> ReedSolomon::decodePGZ(NTL::GF2EX rec,
                                                 int *const resErrs) {
  // maximum correctable errors
  unsigned t = (N - K + 1) / 2;
  NTL::Mat<NTL::GF2E> M;

  std::vector<NTL::GF2E> syndromes(N - K);
  bool allZero = true;
  std::transform(primPow.begin() + primRootOffset,
                 primPow.begin() + primRootOffset + (N - K), syndromes.begin(),
                 [&](const NTL::GF2E &alpha) {
                   NTL::GF2E s = NTL::eval(rec, alpha);
                   allZero = allZero && NTL::IsZero(s);
                   return s;
                 });

  // no nonzero syndromes, r is a valid code
  if (allZero) {
    if (resErrs != nullptr)
      *resErrs = 0;
    return rec;
  }

  NTL::Vec<NTL::GF2E> s;

  for (int v = t; v > 0; --v) {
    NTL::Vec<NTL::GF2E> l;
    // attempting to solve with one less error
    M.SetDims(v, v);
    s.SetLength(v);
    l.SetLength(v);

    // initialize M with syndromes, equation (5.7) from Martini's paper
    // optimize maybe
    // NTL does NOT copy the old elements when doing a resize w different
    // column ct
    for (int i = 0; i < v; ++i) {
      for (int j = 0; j < v; ++j) {
        M[i][j] = syndromes[i + j];
      }
    }

    for (int i = 0; i < v; ++i) {
      s[i] = syndromes[v + i];
    }

    // coefficients of locator polynomial
    NTL::GF2E det;
    NTL::solve(det, M, l, s);

    // system can't be solved
    if (NTL::IsZero(det))
      continue;

    l.SetLength(l.length() + 1);
    l[l.length() - 1] = NTL::GF2E(1);
    std::vector<int> locs = findRootPows(l);

    // failure, couldn't find correct number of error locations
    if (locs.size() < v)
      continue;

    auto vals = solveErrorVals(syndromes, locs);
    if (!vals)
      continue;

    for (int i = 0; i < v; ++i) {
      // final values of added column are exactly the error values
      NTL::SetCoeff(rec, locs[i], NTL::coeff(rec, locs[i]) - vals.value()[i]);
      // NTL::SetCoeff(r, locs[i], NTL::coeff(r, locs[i]) - l[i]);
    }

    if (resErrs != nullptr)
      *resErrs = v;
    return rec;
  }

  // indeterminate number of errors > t, cannot correct
  if (resErrs != nullptr)
    *resErrs = -1;
  return std::nullopt;
}

std::optional<NTL::GF2EX> ReedSolomon::decodeBM(NTL::GF2EX rec,
                                                int *const resErrs) {
  // maximum correctable errors
  unsigned t = (N - K + 1) / 2;
  std::vector<NTL::GF2E> syndromes(N - K);

  for (int i = 0; i < N - K; ++i) {
    syndromes[i] = NTL::eval(rec, primPow[primRootOffset + i]);
  }

  NTL::GF2EX loc_m(1), loc_r(1);
  int m = -1;
  int len_m = 0, len_r = 0;
  NTL::GF2E delta_m(1);
  for (int r = 0; r < (int)(N - K); ++r) {
    NTL::GF2E delta_r(0);
    for (int h = 0; h <= len_r; ++h) {
      delta_r += NTL::coeff(loc_r, h) * syndromes[r - h];
    }

    if (NTL::IsZero(delta_r)) {
      // no change, current solution still still solves for (r + 1)th syndrome
      continue;
    }

    NTL::GF2EX locTmp = loc_r - ((delta_r * NTL::inv(delta_m)) *
                                 (NTL::LeftShift(loc_m, r - m)));
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
  if (NTL::deg(loc_r) > t) {
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

  {
    NTL::GF2EX e;

    for (int i = 0; i < locs.size(); ++i) {
      NTL::SetCoeff(e, locs[i], vals.value()[i]);
    }

    rec -= e;
  }

  // final sanity check, are syndromes zero now?
  for (int i = 0; i < N - K; ++i) {
    if (!NTL::IsZero(NTL::eval(rec, primPow[primRootOffset + i]))) {
      if (resErrs != nullptr)
        *resErrs = -1;
      return std::nullopt;
    }
  }

  if (resErrs != nullptr)
    *resErrs = locs.size();
  return rec;
}

std::vector<int> ReedSolomon::findRootPows(NTL::Vec<NTL::GF2E> v) {
  std::vector<int> res;
  for (int i = 0; i < N; ++i) {
    if (NTL::IsZero(std::accumulate(v.begin(), v.end(), NTL::GF2E::zero()))) {
      res.push_back(i);
    }

    for (int j = 0; j < v.length(); ++j) {
      v[j] *= primPow[j];
    }
  }
  return res;
}

std::vector<int> ReedSolomon::findRootPows(const NTL::GF2EX &loc) {
  NTL::Vec<NTL::GF2E> v = NTL::VectorCopy(NTL::reverse(loc), NTL::deg(loc) + 1);
  std::vector<int> res;
  for (int i = 0; i < N; ++i) {
    if (NTL::IsZero(std::accumulate(v.begin(), v.end(), NTL::GF2E::zero()))) {
      res.push_back(i);
    }

    for (int j = 0; j < v.length(); ++j) {
      v[j] *= primPow[j];
    }
  }
  return res;
}
std::optional<std::vector<NTL::GF2E>>
ReedSolomon::solveErrorVals(const std::vector<NTL::GF2E> &syndromes,
                            const std::vector<int> &errorLocs) {
  assert(syndromes.size() == N - K);
  int v = errorLocs.size();
  NTL::Mat<NTL::GF2E> M(NTL::INIT_SIZE, N - K, v + 1);

  // initialize main M
  for (int i = 0; i < v; ++i) {
    // IMPORTANT: syndrome 0 is defined wrt the used offset!!
    M[0][i] = NTL::power(primPow[errorLocs[i]], primRootOffset);
  }
  for (int i = 1; i < N - K; ++i) {
    for (int j = 0; j < v; ++j) {
      // nasty nasty
      M[i][j] = M[i - 1][j] * primPow[errorLocs[j]];
    }
  }
  for (int i = 0; i < N - K; ++i) {
    M[i][v] = syndromes[i];
  }
  // reduce system
  NTL::gauss(M, v);
  // inconsistent -> no solution
  // dependent -> no unique solution
  if (!checkSystem(M, v))
    return std::nullopt;

  reduce(M, v);

  std::vector<NTL::GF2E> res;
  for (int i = 0; i < v; ++i) {
    res.push_back(M[i][v]);
  }
  return res;
}

std::optional<std::vector<NTL::GF2E>>
ReedSolomon::solveErrorValsForney(const NTL::GF2EX &locator,
                                  const std::vector<NTL::GF2E> &syndromes,
                                  const std::vector<int> &errorLocs) {
  NTL::GF2EX O;

  // initialize with syndrome polynomial
  for (int i = 0; i < syndromes.size(); ++i) {
    NTL::SetCoeff(O, i, syndromes[i]);
  }

  O = NTL::trunc(O * locator, N - K);

  NTL::GF2EX locatorDiff = NTL::diff(locator);

  std::vector<NTL::GF2E> vals;
  vals.resize(errorLocs.size());
  for (int i = 0; i < errorLocs.size(); ++i) {
    NTL::GF2E loc = primPow[errorLocs[i]];
    NTL::GF2E locInv = NTL::inv(loc);

    // Gill, John EE387, source from wikipedia. for non-1 offsets the first
    // term is necessary
    NTL::GF2E l = NTL::eval(locatorDiff, locInv);
    if (NTL::IsZero(l))
      return std::nullopt;

    if (primRootOffset == 0)
      vals[i] = loc * NTL::eval(O, locInv) / l;
    else
      vals[i] = NTL::power(loc, 1 - (int)primRootOffset) *
                     NTL::eval(O, locInv) / l;
  }

  return vals;
}
