#include <stdio.h>
#include <iostream>
#include <complex>

int main(int argc, char ** argv) {
  const std::complex<double> j(0.0, 1.0);
  std::complex<double> z1(-1.0, 1.0);
  std::complex<double> z2(2.0, -2.0);

  printf("Hello World\n");

  std::cout << z1.real() << ", " << z1.imag() << "i" << "\n";
  std::cout << arg(z1) << "\n";
  std::cout << norm(z1) << "\n";
  std::cout << conj(z1) << "\n";
}
