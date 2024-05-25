#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <galois/GaloisField.h>
#include <givaro/gfq.h>
#include <givaro/givmatrix.h>

constexpr int ALPHA = 2;
using SizeStore = int32_t;
using USizeStore = uint32_t;
using GFDom = Givaro::GFqDom<SizeStore>;
using GFPolyDom = Givaro::Poly1Dom<GFDom, Givaro::Dense>;
using GFMatDom = Givaro::MatrixDom<GFDom, Givaro::Dense>;
using GFVecDom = Givaro::VectorDom<GFDom, Givaro::Dense>;

struct ReedSolomon {
  uint Q;
  uint N;
  uint K;
  GFDom gf;
  GFPolyDom gfP;
  GFMatDom gfM;
  GFVecDom gfV;
  // galois field with LUT used for mul and exp
  galois::GaloisField gfLUT;
  uint primPowOffset;

  GFPolyDom::Element genPoly;

  ReedSolomon(uint Q, uint K, const std::vector<USizeStore> &primPoly,
              uint primPowOffset = 0);

  GFPolyDom::Element encode(std::vector<GFDom::Element> &u);

  // returns nullopt if at least N - K errors => error cannot be solved for
  std::optional<GFPolyDom::Element> decodePGZ(const GFPolyDom::Element &r,
                                              size_t const *resErrs = nullptr);

  std::optional<GFPolyDom::Element> decodeBM(const GFPolyDom::Element &r,
                                             size_t const *resErrs = nullptr);
};

