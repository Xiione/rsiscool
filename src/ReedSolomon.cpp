#include <optional>
#include <vector>

#include <galois/GaloisField.h>
#include <galois/GaloisFieldElement.h>
#include <galois/GaloisFieldPolynomial.h>
#include <galois/NumTraits.hpp>

#include "ReedSolomon.hpp"

using namespace galois;

GaloisFieldPolynomial createGenPoly(CodeFormat &fmt) {
  size_t genDeg = fmt.N - fmt.K;
  std::vector<GaloisFieldElement> coeffInit = {
      GaloisFieldElement(&fmt.field, 1)};
  GaloisFieldPolynomial genPoly(&fmt.field, 0, coeffInit.data());

  // (z - a^i)
  std::vector<GaloisFieldElement> minPolyCoeff = {
      GaloisFieldElement(&fmt.field, 0), GaloisFieldElement(&fmt.field, 1)};

  for (size_t i = 0; i < genDeg; ++i) {
    minPolyCoeff[0] = fmt.primEl ^ (fmt.primPowOffset + i);
    // *= (z - a^(i + 1))
    genPoly *= GaloisFieldPolynomial(&fmt.field, 1, minPolyCoeff.data());
  }

  assert(genPoly.deg() == genDeg);
  return genPoly;
}

GaloisFieldPolynomial encode(CodeFormat &fmt,
                             std::vector<GaloisFieldElement> &u) {
  GaloisFieldPolynomial c(&fmt.field, fmt.K - 1, u.data());

  // u(z) * z^(n - k)
  c <<= fmt.N - fmt.K;
  assert(c.deg() == fmt.N - 1);

  // u(z) * z^(n - k) - R_{g(z)}[u(z) * z^(n-k)]
  return c - (c % createGenPoly(fmt));
}

std::optional<GaloisFieldPolynomial> decodePGZ(CodeFormat &fmt,
                                               const GaloisFieldPolynomial &r,
                                               size_t const *res_errs) {
  using namespace Eigen;

  // maximum correctable errors
  size_t t = (fmt.N - fmt.K + 1) / 2;

  Matrix<GaloisFieldElement, Dynamic, Dynamic> M(t - 1, t - 1);
  // initialize M with syndromes, equation (5.7) from Martini's paper
  // optimize maybe
  for (int i = 0; i < t - 1; ++i) {
    for (int j = 0; j < t - 1; ++j) {
      M(i, j) = r(fmt.primEl ^ (i + j));
    }
  }

  Matrix<GaloisFieldElement, Dynamic, 1> S(t - 1);
  for (int i = 0; i < t - 1; ++i) {
    S(i) = r(fmt.primEl ^ (t - 1 + i));
  }

  for (int v = t - 1; v > 0; --v) {
    // attempting to solve with one less error
    M.conservativeResize(v, v);
    S.conservativeResize(v);

    // coefficients of locator polynomial
    Matrix<GaloisFieldElement, Dynamic, 1> V = M.colPivHouseholderQr().solve(S);
    bool solExists = (M * V).isApprox(S);

    std::cout << "v=" << v << ": " << 
  }

  return std::nullopt;
}

std::optional<GaloisFieldPolynomial> decodeBM(CodeFormat &fmt,
                                              const GaloisFieldPolynomial &r,
                                              size_t const *res_errs) {
  return std::nullopt;
}
