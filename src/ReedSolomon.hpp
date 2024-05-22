#pragma once

#include <vector>

#include <galois/GaloisField.h>
#include <galois/GaloisFieldElement.h>
#include <galois/GaloisFieldPolynomial.h>

#define ALPHA 2

galois::GaloisFieldPolynomial createGenPoly(galois::GaloisField &gf, size_t deg,
                                            size_t offset);

galois::GaloisFieldPolynomial
encode(galois::GaloisField &gf, std::vector<galois::GaloisFieldElement> &v_u,
       size_t N);

std::vector<galois::GaloisFieldElement>
decodePGZ(galois::GaloisField &gf, const galois::GaloisFieldPolynomial &r,
       size_t N, size_t K);

std::vector<galois::GaloisFieldElement>
decodeBM(galois::GaloisField &gf, const galois::GaloisFieldPolynomial &r,
       size_t N, size_t K);
