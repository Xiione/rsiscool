#include <vector>

#include <Eigen/Core>
#include <galois/GaloisField.h>
#include <galois/GaloisFieldElement.h>
#include <galois/GaloisFieldPolynomial.h>
#include <galois/NumTraits.hpp>

#include "ReedSolomon.hpp"

using namespace galois;

GaloisFieldPolynomial createGenPoly(GaloisField &gf, size_t deg,
                                    size_t offset) {
  std::vector<GaloisFieldElement> coeffInit = {GaloisFieldElement(&gf, 1)};
  GaloisFieldPolynomial genPoly(&gf, 0, coeffInit.data());

  // (z - a^i)
  std::vector<GaloisFieldElement> minPolyCoeff = {GaloisFieldElement(&gf, 0),
                                                  GaloisFieldElement(&gf, 1)};

  GaloisFieldElement curPrimElPow(&gf, 1);

  for (size_t i = 0; i < offset; ++i) {
    curPrimElPow *= GaloisFieldElement(&gf, ALPHA);
  }

  for (size_t i = 0; i < deg; ++i) {
    minPolyCoeff[0] = curPrimElPow;
    // *= (z - a^(i + 1))
    genPoly *= GaloisFieldPolynomial(&gf, 1, minPolyCoeff.data());
    curPrimElPow *= GaloisFieldElement(&gf, ALPHA);
  }

  assert(genPoly.deg() == deg);
  return genPoly;
}

GaloisFieldPolynomial encode(GaloisField &gf,
                             std::vector<GaloisFieldElement> &v_u, size_t N) {
  size_t K = v_u.size();
  GaloisFieldPolynomial c(&gf, K - 1, v_u.data());

  // u(z) * z^(n - k)
  c <<= N - K;
  assert(c.deg() == N - 1);

  // u(z) * z^(n - k) - R_{g(z)}[u(z) * z^(n-k)]
  return c - (c % createGenPoly(gf, N - K, 1));
}

std::vector<GaloisFieldElement> decodePGZ(GaloisField &gf,
                                           const GaloisFieldPolynomial &r,
                                           size_t N, size_t K) {
  // stuff
}

std::vector<GaloisFieldElement> decodeBM(GaloisField &gf,
                                           const GaloisFieldPolynomial &r,
                                           size_t N, size_t K) {
  // stuff
}
