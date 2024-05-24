#pragma once

#include <optional>
#include <vector>

#include <galois/GaloisField.h>
#include <galois/GaloisFieldElement.h>
#include <galois/GaloisFieldPolynomial.h>

#define ALPHA 2

struct CodeFormat {
  galois::GaloisField &field;
  size_t N;
  size_t K;
  galois::GaloisFieldElement primEl;
  size_t primPowOffset;

  CodeFormat(galois::GaloisField &field, size_t N, size_t K,
             size_t primPowOffset = 0)
      : field(field), N(N), K(K), primPowOffset(primPowOffset) {
    primEl = galois::GaloisFieldElement(&field, ALPHA);
  }
};

galois::GaloisFieldPolynomial createGenPoly(CodeFormat &fmt);

galois::GaloisFieldPolynomial
encode(CodeFormat &fmt, std::vector<galois::GaloisFieldElement> &u);

// returns nullopt if at least N - K errors => error cannot be solved for
std::optional<galois::GaloisFieldPolynomial>
decodePGZ(CodeFormat &fmt, const galois::GaloisFieldPolynomial &r,
          size_t const *resErrs = nullptr);

std::optional<galois::GaloisFieldPolynomial>
decodeBM(CodeFormat &fmt, const galois::GaloisFieldPolynomial &r,
          size_t const *resErrs = nullptr);
