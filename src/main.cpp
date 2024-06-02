#include <cstdint>
#include <iostream>
#include <vector>

#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <galois/GaloisFieldElement.h>

#include "ReedSolomon.hpp"

int main() {
  std::vector<uint> p = {1, 1, 0, 1};
  galois::GaloisField gf(3, p.data());
  galois::GaloisFieldElement el(&gf, 2);
  auto elp = el ^ 0;
  auto elp2 = el ^ 1;

  const std::vector<uint> symbolPrimPoly = {1, 1, 0, 1};
  NTL::GF2X primPoly = vecToGF2X(symbolPrimPoly);
  std::cout << "Using primitive polynomial: " << primPoly << std::endl;
  NTL::GF2E::init(primPoly);

  ReedSolomon rs(8, 3, 0);

  std::cout << "Enter u:" << std::endl;

  std::vector<NTL::GF2E> u;
  for (uint i = 0; i < rs.K; ++i) {
    int x;
    std::cin >> x;
    u.push_back(NTL::to_GF2E(intToGF2X(x)));
    std::cout << u.back() << std::endl;
  }

  NTL::GF2EX c = rs.encode(u);

  std::cout << "Encoded c:" << std::endl;
  for (uint i = 0; i <= NTL::deg(c); ++i)
    std::cout << GF2EtoInt(NTL::coeff(c, i)) << ' ';
  std::cout << std::endl;

  std::cout << "Enter received c (r):" << std::endl;
  NTL::GF2EX r;
  for (uint i = 0; i < rs.N; ++i) {
    int x;
    std::cin >> x;
    NTL::SetCoeff(r, i, intToGF2E(x));
  }

  std::cout << "r: " << std::endl;
  for (uint i = 0; i <= NTL::deg(r); ++i)
    std::cout << GF2EtoInt(NTL::coeff(r, i)) << ' ';
  std::cout << std::endl;

  int errsPGZ, errsBM;
  auto start = std::chrono::high_resolution_clock::now();
  auto resPGZ = rs.decodePGZ(r, &errsPGZ);
  auto between = std::chrono::high_resolution_clock::now();
  auto resBM = rs.decodeBM(r, &errsBM);
  auto end = std::chrono::high_resolution_clock::now();

  std::cout << "PGZ Decoder took "
            << std::chrono::duration_cast<std::chrono::nanoseconds>(between -
                                                                    start)
                   .count()
            << " nanoseconds:" << std::endl;

  if (resPGZ) {
    std::cout << "r with " << errsPGZ << " error(s) corrected: " << std::endl;
    for (uint i = 0; i <= NTL::deg(*resPGZ); ++i)
      std::cout << GF2EtoInt(NTL::coeff(*resPGZ, i)) << ' ';
    std::cout << std::endl;
  } else
    std::cout << "PGZ decoding failed" << std::endl;
  std::cout << "BM Decoder took "
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                    between)
                   .count()
            << " nanoseconds:" << std::endl;
  if (resBM) {
    std::cout << "r with " << errsBM << " error(s) corrected: " << std::endl;
    for (uint i = 0; i <= NTL::deg(*resBM); ++i)
      std::cout << GF2EtoInt(NTL::coeff(*resBM, i)) << ' ';
    std::cout << std::endl;
  } else
    std::cout << "BM decoding failed" << std::endl;
}
