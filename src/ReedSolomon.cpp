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

ReedSolomon::ReedSolomon(uint Q, uint K, uint offset)
    : Q(Q), N(Q - 1), K(K), primRootOffset(offset) {

  // deduce extension field power
  assert(Q == 1 << (int)log2(Q));

  // init generator polynomial
  uint genDeg = N - K;
  NTL::GF2EX primMinPoly;

  NTL::GF2E a = intToGF2E(ALPHA), x(1);
  assert(NTL::IsOne(x) != 0);
  for (uint i = 0; i < N; ++i) {
    primPow.push_back(x);
    x *= a;
  }

  // (z - a^i)
  NTL::Vec<NTL::GF2E> roots;
  roots.SetLength(genDeg);
  for (size_t i = 0; i < genDeg; ++i) {
    roots[i] = primPow[offset + i];
  }

  genPoly = NTL::BuildFromRoots(roots);
  assert(NTL::deg(genPoly) == genDeg);
}

NTL::GF2EX ReedSolomon::encode(std::vector<NTL::GF2E> &u) {
  NTL::GF2EX c;
  for (uint i = 0; i < u.size(); ++i)
    NTL::SetCoeff(c, i, u[i]);

  // u(z) * z^(n - k)
  c <<= N - K;
  assert(NTL::deg(c) == N - 1);

  // u(z) * z^(n - k) - R_{g(z)}[u(z) * z^(n-k)]
  c -= c % genPoly;
  return c;
}

std::optional<NTL::GF2EX> ReedSolomon::decodePGZ(NTL::GF2EX rec,
                                                 int *const resErrs) {
  // maximum correctable errors
  unsigned t = (N - K + 1) / 2;
  NTL::Mat<NTL::GF2E> M;

  std::vector<NTL::GF2E> syndromes(N - K);

  bool allZero = true;
  for (uint i = 0; i < N - K; ++i) {
    syndromes[i] = NTL::eval(rec, primPow[primRootOffset + i]);
    allZero = allZero && NTL::IsZero(syndromes[i]);
  }

  // no nonzero syndromes, r is a valid code
  if (allZero) {
    if (resErrs != nullptr)
      *resErrs = 0;
    return rec;
  }

  NTL::Vec<NTL::GF2E> s;

  for (uint v = t; v > 0; --v) {
    NTL::Vec<NTL::GF2E> l;
    // attempting to solve with one less error
    M.SetDims(v, v);
    s.SetLength(v);
    l.SetLength(v);

    // initialize M with syndromes, equation (5.7) from Martini's paper
    // optimize maybe
    // NTL does NOT copy the old elements when doing a resize w different
    // column ct
    for (uint i = 0; i < v; ++i) {
      for (uint j = 0; j < v; ++j) {
        M[i][j] = syndromes[i + j];
      }
    }

    for (uint i = 0; i < v; ++i) {
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
    std::vector<uint> locs = findRootPows(l);

    // failure, couldn't find correct number of error locations
    if (locs.size() < v)
      continue;

    auto vals = solveErrorVals(syndromes, locs);
    if (!vals)
      continue;

    for (uint i = 0; i < v; ++i) {
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
  // unsigned t = (N - K + 1) / 2;
  std::vector<NTL::GF2E> syndromes(N - K);

  for (uint i = 0; i < N - K; ++i) {
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

  std::vector<uint> locs =
      findRootPows(NTL::VectorCopy(NTL::reverse(loc_r), NTL::deg(loc_r) + 1));

  auto vals = solveErrorValsForney(loc_r, syndromes, locs);

  if (!vals) {
    if (resErrs != nullptr)
      *resErrs = -1;
    return std::nullopt;
  }

  if (resErrs != nullptr)
    *resErrs = locs.size();

  for (uint i = 0; i < locs.size(); ++i) {
    NTL::SetCoeff(rec, locs[i], NTL::coeff(rec, locs[i]) - vals.value()[i]);
  }

  return rec;
}

std::vector<uint> ReedSolomon::findRootPows(NTL::Vec<NTL::GF2E> v) {
  std::unordered_set<uint> s;
  for (uint i = 0; i < N; ++i) {
    if (NTL::IsZero(std::accumulate(v.begin(), v.end(), NTL::GF2E(0)))) {
      s.insert(i);
    }

    for (uint j = 0; j < v.length(); ++j) {
      v[j] *= primPow[j];
    }
  }
  return std::vector<uint>(s.begin(), s.end());
}

std::optional<std::vector<NTL::GF2E>>
ReedSolomon::solveErrorVals(const std::vector<NTL::GF2E> &syndromes,
                            const std::vector<uint> &locs) {
  assert(syndromes.size() == N - K);
  uint v = locs.size();
  NTL::Mat<NTL::GF2E> M(NTL::INIT_SIZE, N - K, v + 1);

  // initialize main M
  for (uint i = 0; i < v; ++i) {
    // IMPORTANT: syndrome 0 is defined wrt the used offset!!
    M[0][i] = NTL::power(primPow[locs[i]], primRootOffset);
  }
  for (uint i = 1; i < N - K; ++i) {
    for (uint j = 0; j < v; ++j) {
      // nasty nasty
      M[i][j] = M[i - 1][j] * primPow[locs[j]];
    }
  }
  for (uint i = 0; i < N - K; ++i) {
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
  for (uint i = 0; i < v; ++i) {
    res.push_back(M[i][v]);
  }
  return res;
}

std::optional<std::vector<NTL::GF2E>>
ReedSolomon::solveErrorValsForney(const NTL::GF2EX &locator,
                                  const std::vector<NTL::GF2E> &syndromes,
                                  const std::vector<uint> &locs) {
  assert(syndromes.size() == N - K);
  NTL::GF2EX O;

  // initialize with syndrome polynomial
  for (uint i = 0; i < syndromes.size(); ++i) {
    NTL::SetCoeff(O, i, syndromes[i]);
  }

  O = NTL::trunc(O * locator, N - K);

  NTL::GF2EX locatorDiff = NTL::diff(locator);

  std::vector<NTL::GF2E> vals;
  for (uint i = 0; i < locs.size(); ++i) {
    NTL::GF2E loc = primPow[locs[i]];
    NTL::GF2E locInv = NTL::inv(loc);

    // Gill, John EE387, source from wikipedia. for non-1 offsets the first
    // term is necessary
    vals.push_back(NTL::power(loc, 1 - (int)primRootOffset) * NTL::eval(O, locInv) /
                   NTL::eval(locatorDiff, locInv));
  }

  return vals;
}

// random
NTL::GF2X intToGF2X(uint x) {
  size_t numBytes = (sizeof(uint) * 8 - __builtin_clz(x) + 7) / 8;
  return NTL::GF2XFromBytes((uint8_t *)&x, numBytes);
}

NTL::GF2E intToGF2E(uint x) { return NTL::to_GF2E(intToGF2X(x)); }

NTL::GF2X vecToGF2X(const std::vector<uint> v) {
  std::vector<uint8_t> res;
  size_t numBytes = (v.size() + 7) / 8;

  for (size_t i = 0; i < numBytes; ++i) {
    uint8_t byte = 0;
    for (size_t j = 0; j < 8; ++j) {
      size_t index = i * 8 + j;
      if (index < v.size()) {
        if (v[index] != 0) {
          byte |= (1 << j);
        }
      } else
        break;
    }
    res.push_back(byte);
  }
  return NTL::GF2XFromBytes(res.data(), numBytes);
}

uint GF2EtoInt(const NTL::GF2E &x) {
  NTL::GF2X f = x._GF2E__rep;
  uint res = 0;
  for (uint i = 0; i <= NTL::deg(f); ++i) {
    if (!NTL::IsZero(NTL::coeff(f, i)))
      res |= 1 << i;
  }
  return res;
}
