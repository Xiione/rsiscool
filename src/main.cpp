#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <galois/GaloisField.h>
#include <galois/GaloisFieldElement.h>
#include <galois/GaloisFieldPolynomial.h>
#include <galois/NumTraits.hpp>

#include "ReedSolomon.hpp"

constexpr int Q = 8;
constexpr int N = Q - 1;
constexpr int K = 3;

// const std::vector<unsigned int> BIGSYMBOL_PRIMPOLY = {1, 1, 1, 1, 0, 0, 0,
// 1};
const std::vector<unsigned int> SYMBOL_PRIMPOLY = {1, 1, 0, 1};
galois::GaloisField SYMBOL_GF(32 - __builtin_clz(Q) - 1,
                              SYMBOL_PRIMPOLY.data());
const galois::GaloisFieldElement SYMBOL_PRIMEL(&SYMBOL_GF, ALPHA);

template <typename Func>
auto measure_time(Func func, const std::vector<int> &data) {
  auto start = std::chrono::high_resolution_clock::now();
  func(data);
  auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
      .count();
}

int main() {
  std::cout << "Enter u:" << std::endl;
  std::vector<galois::GaloisFieldElement> u;
  for (int i = 0; i < K; ++i) {
    int x;
    std::cin >> x;
    u.emplace_back(&SYMBOL_GF, x);
  }

  galois::GaloisFieldPolynomial c = encode(SYMBOL_GF, u, N);

  std::cout << "Encoded c:" << std::endl;
  for (int i = 0; i < N; ++i) {
    std::cout << c[i].poly() << ' ';
  }
  std::cout << std::endl;

  std::cout << "Enter received c:" << std::endl;
  std::vector<galois::GaloisFieldElement> vr;
  for (int i = 0; i < K; ++i) {
    int x;
    std::cin >> x;
    u.emplace_back(&SYMBOL_GF, x);
  }

  galois::GaloisFieldPolynomial r(&SYMBOL_GF, N, vr.data());

  auto start = std::chrono::high_resolution_clock::now();
  auto resPGZ = decodePGZ(SYMBOL_GF, r, N, K);
  auto between = std::chrono::high_resolution_clock::now();
  auto resBM = decodeBM(SYMBOL_GF, r, N, K);
  auto end = std::chrono::high_resolution_clock::now();

  std::cout << "PGZ Decoder took "
            << std::chrono::duration_cast<std::chrono::nanoseconds>(between -
                                                                    start)
                   .count() << " nanoseconds:" << std::endl;
  for(auto &gfe : resPGZ)
    std::cout << gfe.poly() << ' ';
  std::cout << std::endl;
  std::cout << "BM Decoder took "
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                    between)
                   .count() << " nanoseconds:" << std::endl;
  for(auto &gfe : resBM)
    std::cout << gfe.poly() << ' ';
  std::cout << std::endl;
}
