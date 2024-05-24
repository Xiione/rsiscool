#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <galois/GaloisField.h>
#include <galois/GaloisFieldElement.h>
#include <galois/GaloisFieldPolynomial.h>
#include <galois/NumTraits.hpp>

#include "ReedSolomon.hpp"

template <typename Func>
auto measure_time(Func func, const std::vector<int> &data) {
  auto start = std::chrono::high_resolution_clock::now();
  func(data);
  auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
      .count();
}

int main() {
  int Q = 8;
  const std::vector<unsigned int> symbolPrimPoly = {1, 1, 0, 1};
  galois::GaloisField symbolGf(32 - __builtin_clz(Q) - 1,
                                symbolPrimPoly.data());
  CodeFormat fmt(symbolGf, Q - 1, 3, 1);
  std::cout << "Enter u:" << std::endl;


  std::vector<galois::GaloisFieldElement> u;
  for (int i = 0; i < fmt.K; ++i) {
    int x;
    std::cin >> x;
    u.emplace_back(&symbolGf, x);
  }

  galois::GaloisFieldPolynomial c = encode(fmt, u);

  std::cout << "Encoded c:" << std::endl;
  for (int i = 0; i < fmt.N; ++i) {
    std::cout << c[i].poly() << ' ';
  }
  std::cout << std::endl;

  std::cout << "Enter received c:" << std::endl;
  std::vector<galois::GaloisFieldElement> vr;
  for (int i = 0; i < fmt.N; ++i) {
    int x;
    std::cin >> x;
    u.emplace_back(&symbolGf, x);
  }

  galois::GaloisFieldPolynomial r(&symbolGf, fmt.N, vr.data());

  auto start = std::chrono::high_resolution_clock::now();
  auto resPGZ = decodePGZ(fmt, r);
  auto between = std::chrono::high_resolution_clock::now();
  auto resBM = decodeBM(fmt, r);
  auto end = std::chrono::high_resolution_clock::now();

  std::cout << "PGZ Decoder took "
            << std::chrono::duration_cast<std::chrono::nanoseconds>(between -
                                                                    start)
                   .count()
            << " nanoseconds:" << std::endl;

  if (resPGZ) {
    for (size_t i = 0; i <= resPGZ->deg(); ++i)
      std::cout << resPGZ.value()[i] << ' ';
    std::cout << std::endl;
  } else
    std::cout << "PGZ decoding failed" << std::endl;
  std::cout << "BM Decoder took "
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                    between)
                   .count()
            << " nanoseconds:" << std::endl;
  if (resBM) {
    for (size_t i = 0; i <= resBM->deg(); ++i)
      std::cout << resBM.value()[i] << ' ';
    std::cout << std::endl;
  } else
    std::cout << "BM decoding failed" << std::endl;
}
