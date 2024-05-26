#include <cassert>
#include <optional>
#include <vector>

#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2E.h>
#include <galois/GaloisFieldElement.h>

#include "ReedSolomon.hpp"

ReedSolomon::ReedSolomon(uint Q, uint K, const std::vector<uint> &primPoly,
                         uint primPowOffset)
    : Q(Q), N(Q - 1), K(K), gfLUT(log2(Q), primPoly.data()),
      primPowOffset(primPowOffset) {

  // deduce extension field power
  assert(Q == 1 << (int)log2(Q));

  // init generator polynomial
  uint genDeg = N - K;
  NTL::GF2EX primMinPoly;

  NTL::SetCoeff(genPoly, 0, 1);
  NTL::SetCoeff(primMinPoly, 0, 1);
  NTL::SetCoeff(primMinPoly, 1, 1);
  // (z - a^i)
  for (size_t i = 0; i < genDeg; ++i) {
    galois::GaloisFieldElement el(&gfLUT, ALPHA);
    auto elp = el ^ (0);
    uint rep = (galois::GaloisFieldElement(&gfLUT, ALPHA) ^ (primPowOffset + i))
                   .poly();
    NTL::GF2E coeff = intToGF2E(rep);
    NTL::SetCoeff(primMinPoly, 0, coeff);
    // dom_power(primMinPoly[0], (GFDom::Element)gf.generator(), primPowOffset +
    // i, gf);
    genPoly *= primMinPoly;
  }
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

std::optional<NTL::GF2EX> ReedSolomon::decodePGZ(const NTL::GF2EX &r,
                                                 size_t const *resErrs) {
  // maximum correctable errors
  unsigned t = (N - K + 1) / 2;
  NTL::Mat<NTL::GF2E> M;
  M.SetDims(t - 1, t - 1);

  // initialize M with syndromes, equation (5.7) from Martini's paper
  // optimize maybe
  std::vector<NTL::GF2E> syndromes;
  for (uint i = 0; i < 2 * t - 1; ++i) {
    NTL::GF2E val =
        intToGF2E((galois::GaloisFieldElement(&gfLUT, ALPHA) ^ i).poly());
    syndromes.push_back(val);
  }

  for (uint i = 0; i < t - 1; ++i) {
    for (uint j = 0; j < t - 1; ++j) {
      // M(i, j) = r(fmt.primEl ^ (i + j));
      M(i, j) = NTL::eval(r, syndromes[i + j]);
    }
  }

  NTL::Vec<NTL::GF2E> S;
  S.SetLength(t - 1);

  for (uint i = 0; i < t - 1; ++i) {
    S(i) = NTL::eval(r, syndromes[t - 1 + i]);
  }

  for (int v = t - 1; v > 0; --v) {
    // attempting to solve with one less error
    M.SetDims(v, v);
    S.SetLength(v);

    // coefficients of locator polynomial
  }

  return std::nullopt;
}

std::optional<NTL::GF2EX> ReedSolomon::decodeBM(const NTL::GF2EX &r,
                                                size_t const *resErrs) {
  return std::nullopt;
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
    for (size_t bit = 0; bit < 8; ++bit) {
      size_t index = i * 8 + bit;
      if (index < v.size()) {
        if (v[index] != 0) {
          byte |= (1 << bit);
        }
      } else
        break;
    }
    res.push_back(byte);
  }
  return NTL::GF2XFromBytes(res.data(), numBytes);
}
