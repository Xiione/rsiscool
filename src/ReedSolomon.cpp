#include <optional>
#include <vector>

#include <galois/GaloisFieldElement.h>
#include <givaro/gfq.h>
#include <givaro/givmatrix.h>

#include "ReedSolomon.hpp"

using namespace Givaro;
using namespace galois;

ReedSolomon::ReedSolomon(uint Q, uint K,
                         const std::vector<USizeStore> &primPoly,
                         uint primPowOffset)
    : Q(Q), N(Q - 1), K(K), primPowOffset(primPowOffset) {

  // deduce extension field power
  int e = log2(Q);
  assert(Q == 1 << e);
  gf = GFDom(2, e, primPoly);
  gfP = GFPolyDom(gf, 'z');
  gfM = GFMatDom(gf);
  gfV = GFVecDom(gf);
  gfLUT = galois::GaloisField(e, primPoly.data());

  // init generator polynomial
  uint genDeg = N - K;
  GFPolyDom::Element primMinPoly;
  gfP.init(genPoly, {gf.one});
  gfP.init(primMinPoly, {gf.one, gf.one});
  // (z - a^i)
  for (size_t i = 0; i < genDeg; ++i) {
    primMinPoly[0] =
        (GaloisFieldElement(&gfLUT, gf.generator()) ^ (primPowOffset + i))
            .poly();
    // dom_power(primMinPoly[0], (GFDom::Element)gf.generator(), primPowOffset +
    // i, gf);
    gfP.mulin(genPoly, primMinPoly);
  }
  assert(gfP.degree(genPoly) == genDeg);
}

GFPolyDom::Element ReedSolomon::encode(std::vector<GFDom::Element> &u) {
  GFPolyDom::Element c;
  gfP.init(c, u);

  // u(z) * z^(n - k)
  gfP.shiftin(c, N - K);

  assert(gfP.degree(c) == N - 1);

  // u(z) * z^(n - k) - R_{g(z)}[u(z) * z^(n-k)]
  GFPolyDom::Element tmp;
  gfP.mod(tmp, c, genPoly);
  gfP.subin(c, tmp);

  return c;
}

std::optional<GFPolyDom::Element>
ReedSolomon::decodePGZ(const GFPolyDom::Element &r, size_t const *resErrs) {
  // maximum correctable errors
  unsigned t = (N - K + 1) / 2;
  GFMatDom::Element M;
  gfM.init(M, t - 1, t - 1);

  // initialize M with syndromes, equation (5.7) from Martini's paper
  // optimize maybe
  for (uint i = 0; i < t - 1; ++i) {
    for (uint j = 0; j < t - 1; ++j) {
      // M(i, j) = r(fmt.primEl ^ (i + j));
      gfP.eval(M(i, j), r,
               (GaloisFieldElement(&gfLUT, gf.generator()) ^ (i + j)).poly());
    }
  }

  GFVecDom::Element S;
  gfV.init(S, t - 1);

  for (uint i = 0; i < t - 1; ++i) {
    gfP.eval(S[i], r,
             (GaloisFieldElement(&gfLUT, gf.generator()) ^ (t - 1 + i)).poly());
  }

  for (int v = t - 1; v > 0; --v) {
    // attempting to solve with one less error
    M.resize(v, v);
    S.resize(v);

    // coefficients of locator polynomial
    
  }

  return std::nullopt;
}

std::optional<GFPolyDom::Element>
ReedSolomon::decodeBM(const GFPolyDom::Element &r, size_t const *resErrs) {
  return std::nullopt;
}
