#pragma once

#include "GaloisFieldElement.h"

#include <Eigen/Core>

namespace Eigen {

template <> struct NumTraits<galois::GaloisFieldElement> {
  typedef galois::GaloisFieldElement Real;
  typedef galois::GaloisFieldElement NonInteger;
  typedef galois::GaloisFieldElement Nested;

  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 0,
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 1, // just an xor
    MulCost = 1, // uses LUT by default
  };
};

} // namespace Eigen

namespace galois {

inline const GaloisFieldElement &abs(const GaloisFieldElement &gfe) {
  return gfe;
}
inline GaloisFieldElement abs2(const GaloisFieldElement &gfe) {
  return gfe * gfe;
}
inline GaloisFieldElement pow(const GaloisFieldElement &gfe, const int &n) {
  return gfe ^ n;
}

} // namespace galois
