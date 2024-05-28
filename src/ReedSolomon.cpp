#include <cassert>
#include <numeric>
#include <optional>
#include <vector>

#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2E.h>
#include <galois/GaloisFieldElement.h>

#include "ReedSolomon.hpp"

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

std::optional<NTL::GF2EX> ReedSolomon::decodePGZ(NTL::GF2EX r,
                                                 uint *const resErrs) {
  // maximum correctable errors
  unsigned t = (N - K + 1) / 2;
  NTL::Mat<NTL::GF2E> M;

  NTL::Vec<NTL::GF2E> syndromes;
  syndromes.SetLength(N - K);

  bool allZero = true;
  for (uint i = 0; i < N - K; ++i) {
    syndromes[i] = NTL::eval(r, primPow[i + primRootOffset]);
    allZero = allZero && NTL::IsZero(syndromes[i]);
  }

  // no nonzero syndromes, r is valid
  if (allZero) {
    if (resErrs != nullptr)
      *resErrs = 0;
    return r;
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
    std::vector<uint> locs = findRoots(l);

    // solve for error values using syndromes and newly found error locs
    // hijack mat and vec variables cuz we're done with them
    M.SetDims(N - K, v + 1);
    l.SetLength(v);

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

    // initialize augmented matrix column with syndromes
    for (uint i = 0; i < N - K; ++i) {
      M[i][v] = syndromes[i];
    }

    // reduce system
    NTL::gauss(M, v);
    for (uint i = 0; i < v; ++i) 
      assert(NTL::IsZero(M[v][i]));
    reduce(M, v);

    // move syndromes back after row operations from gaussian elimination
    // syndromes.SetLength(v);
    // for (uint i = 0; i < v; ++i) {
    //   syndromes[i] = M[i][v];
    // }

    // NTL::Mat<NTL::GF2E> X(NTL::INIT_SIZE, v, v);
    // copy contents over for solving
    // for (uint i = 0; i < v; ++i) {
    //   for (uint j = 0; j < v; ++j) {
    //     X[i][j] = M[i][j];
    //   }
    // }

    // NTL::solve(det, X, l, syndromes);
    // assert(!NTL::IsZero(det));

    for (uint i = 0; i < v; ++i) {
      NTL::SetCoeff(r, locs[i], NTL::coeff(r, locs[i]) - M[i][v]);
      // NTL::SetCoeff(r, locs[i], NTL::coeff(r, locs[i]) - l[i]);
    }

    if (resErrs != nullptr)
      *resErrs = v;
    return r;
  }

  return std::nullopt;
}

std::optional<NTL::GF2EX> ReedSolomon::decodeBM(NTL::GF2EX r,
                                                uint *const resErrs) {
  return std::nullopt;
}

std::vector<uint> ReedSolomon::findRoots(NTL::Vec<NTL::GF2E> v) {
  std::vector<uint> res;
  for (uint i = 0; i < N; ++i) {
    if (NTL::IsZero(std::accumulate(v.begin(), v.end(), NTL::GF2E(0)))) {
      res.push_back(i);
    }

    for (uint j = 0; j < v.length(); ++j) {
      v[j] *= primPow[j];
    }
  }
  return res;
}

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

void reduce(NTL::Mat<NTL::GF2E> &A, int rows) {
  for (int i = rows - 1; i >= 0; --i) {
    A[i] *= NTL::inv(A[i][i]);

    for (int j = 0; j < i; ++j) {
      A[j] -= A[j][i] * A[i];
    }
  }


}
