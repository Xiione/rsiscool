#include <iostream>
#include <vector>

#include <galois/GaloisField.h>
#include <galois/GaloisFieldElement.h>
#include <galois/GaloisFieldPolynomial.h>

#include <galois/NumTraits.hpp>

#include <Eigen/Core>

const int Q = 8;
const int N = Q - 1;
const int K = 3;

// const std::vector<unsigned int> SYMBOL_PRIMPOLY = {1, 1, 1, 1, 0, 0, 0, 1};
const std::vector<unsigned int> SYMBOL_PRIMPOLY = {1, 1, 0, 1};
galois::GaloisField SYMBOL_GF(32 - __builtin_clz(Q) - 1,
                              SYMBOL_PRIMPOLY.data());
const galois::GaloisFieldElement SYMBOL_PRIMEL(&SYMBOL_GF, 2);

galois::GaloisFieldElement GFE(galois::GaloisField *gf,
                               galois::GFSymbol value) {
  return galois::GaloisFieldElement(gf, value);
}

int main() {
  std::vector<galois::GaloisFieldElement> code_genpoly_coeff_init = {
      GFE(&SYMBOL_GF, 1)};

  galois::GaloisFieldPolynomial code_genpoly(&SYMBOL_GF, 0,
                                             code_genpoly_coeff_init.data());

  // create code generator polynomial
  {
    std::vector<galois::GaloisFieldElement> minpoly_coeff = {
        GFE(&SYMBOL_GF, 0), GFE(&SYMBOL_GF, 1)};

    galois::GaloisFieldElement cur_primel_pow(&SYMBOL_GF, 1);
    for (int i = 0; i < N - K; ++i) {
      // a^i * a
      cur_primel_pow *= SYMBOL_PRIMEL;
      // std::cout << "Primitive element power " << i + 1 << ": " <<
      // cur_primel_pow
      //           << std::endl;
      minpoly_coeff[0] = cur_primel_pow;
      // std::cout << "Additive inverse: " << GFE(&SYMBOL_GF, 0) -
      // cur_primel_pow
      //           << std::endl;

      // *= (z - a^(i + 1))
      code_genpoly *=
          galois::GaloisFieldPolynomial(&SYMBOL_GF, 1, minpoly_coeff.data());
    }
  }

  std::cout << "Enter u:" << std::endl;

  std::vector<galois::GaloisFieldElement> u;
  for (int i = 0; i < K; ++i) {
    int x;
    std::cin >> x;
    u.emplace_back(&SYMBOL_GF, x);
  }

  galois::GaloisFieldPolynomial c(&SYMBOL_GF, K - 1, u.data());
  // u(z) * z^(n - k)
  c <<= N - K;

  // u(z) * z^(n - k) - R_{g(z)}[u(z) * z^(n-k)]
  c -= (c % code_genpoly);

  std::cout << "Encoded c:" << std::endl;
  for (int i = 0; i < N; ++i) {
    std::cout << c[i] << ' ';
  }
  std::cout << std::endl;

  std::cout << "Enter received c:" << std::endl;

  // std::vector<galois::GaloisFieldElement> c;
  // for (int i = 0; i < K; ++i) {
  //   int x;
  //   std::cin >> x;
  //   u.emplace_back(&SYMBOL_GF, x);
  // }

  Eigen::Matrix<galois::GaloisFieldElement, Eigen::Dynamic, Eigen::Dynamic> m;
}
