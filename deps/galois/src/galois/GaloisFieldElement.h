/*
 *********************************************************************
 *                                                                   *
 *        Galois Field Arithmetic Library (version 0.0.1)            *
 *                                                                   *
 * Class: Galois Field Element                                       *
 * Version: 0.0.1                                                    *
 * Author: Arash Partow - 2000                                       *
 * URL: http://www.partow.net/projects/galois/index.html             *
 *                                                                   *
 * Copyright Notice:                                                 *
 * Free use of this library is permitted under the guidelines and    *
 * in accordance with the most current version of the Common Public  *
 * License.                                                          *
 * http://www.opensource.org/licenses/cpl.php                        *
 *                                                                   *
 *********************************************************************
 */

#ifndef INCLUDE_GALOISFIELDELEMENT_H
#define INCLUDE_GALOISFIELDELEMENT_H

#include "GaloisField.h"
#include "GlobalGaloisField.h"
#include <iostream>
#include <vector>

namespace galois {

class GaloisFieldElement {

public:
  GaloisFieldElement(GFSymbol v = -1);
  GaloisFieldElement(const GaloisFieldElement &gfe);
  ~GaloisFieldElement() {}

  inline GaloisFieldElement &operator=(const GaloisFieldElement &gfe) {
    if (this == &gfe)
      return *this;

    poly_value = gfe.poly_value;

    return *this;
  }

  inline GaloisFieldElement &operator=(const GFSymbol &v) {
    poly_value = v & GF.size();
    return *this;
  }

  inline GaloisFieldElement &operator+=(const GaloisFieldElement &gfe) {
    poly_value ^= gfe.poly_value;
    return *this;
  }

  inline GaloisFieldElement &operator+=(const GFSymbol &v) {
    poly_value ^= v;
    return *this;
  }

  inline GaloisFieldElement &operator-=(const GaloisFieldElement &gfe) {
    *this += gfe;
    return *this;
  }

  inline GaloisFieldElement &operator-=(const GFSymbol &v) {
    *this += v;
    return *this;
  }

  inline GaloisFieldElement &operator*=(const GaloisFieldElement &gfe) {
    poly_value = GF.mul(poly_value, gfe.poly_value);
    return *this;
  }

  inline GaloisFieldElement &operator*=(const GFSymbol &v) {
    poly_value = GF.mul(poly_value, v);
    return *this;
  }

  inline GaloisFieldElement &operator/=(const GaloisFieldElement &gfe) {
    poly_value = GF.div(poly_value, gfe.poly_value);
    return *this;
  }

  inline GaloisFieldElement &operator/=(const GFSymbol &v) {
    poly_value = GF.div(poly_value, v);
    return *this;
  }

  inline GaloisFieldElement &operator^=(const int &n) {
    poly_value = GF.exp(poly_value, n);
    return *this;
  }

  inline bool operator==(const GaloisFieldElement &gfe) const {
    return (poly_value == gfe.poly_value);
  }

  inline bool operator==(const GFSymbol &v) const { return (poly_value == v); }

  inline bool operator!=(const GaloisFieldElement &gfe) const {
    return (poly_value != gfe.poly_value);
  }

  inline bool operator!=(const GFSymbol &v) const { return (poly_value != v); }

  inline bool operator<(const GaloisFieldElement &gfe) {
    return (poly_value < gfe.poly_value);
  }

  inline bool operator<(const GFSymbol &v) { return (poly_value < v); }

  inline bool operator>(const GaloisFieldElement &gfe) {
    return (poly_value > gfe.poly_value);
  }

  inline bool operator>(const GFSymbol &v) { return (poly_value > v); }

  inline GFSymbol index() const { return GF.index(poly_value); }

  inline GFSymbol poly() const { return poly_value; }

  inline GFSymbol inverse() const { return GF.inverse(poly_value); }

  friend std::ostream &operator<<(std::ostream &os,
                                  const GaloisFieldElement &gfe);

private:
  GFSymbol poly_value;
};

GaloisFieldElement operator+(const GaloisFieldElement &a,
                             const GaloisFieldElement &b);
GaloisFieldElement operator-(const GaloisFieldElement &a,
                             const GaloisFieldElement &b);
GaloisFieldElement operator*(const GaloisFieldElement &a,
                             const GaloisFieldElement &b);
GaloisFieldElement operator*(const GaloisFieldElement &a, const GFSymbol &b);
GaloisFieldElement operator*(const GFSymbol &a, const GaloisFieldElement &b);
GaloisFieldElement operator/(const GaloisFieldElement &a,
                             const GaloisFieldElement &b);
GaloisFieldElement operator^(const GaloisFieldElement &a, const int &b);

} // namespace galois

#endif
