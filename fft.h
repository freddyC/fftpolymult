#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex>
#include <cmath>
#include <utility>

namespace fft
{
  const double PI = 3.14159265359;

  std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > getOmegasForSize (int n) {
    std::vector<std::complex<double> > omega(n);
    std::vector<std::complex<double> > invOmega(n);

    for (int i = 0; i < n; ++i) {
      omega[i] = std::complex<double> (std::cos((2 * PI * i)/n), std::sin((2 * PI * i)/n));
      invOmega[i] = std::complex<double> (std::cos((2 * PI * i)/n), -1 * std::sin((2 * PI * i)/n));
    }
    return std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > (omega, invOmega);
  }

  std::vector<std::complex<double> > fftRecursive(std::vector<std::complex<double> > &poly, int power, bool useOmega, std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > &omegas) {
    // BASE CASE
    if (poly.size() == 1) {
      return poly;
    }

    // Divide into two parts
    int n = poly.size();
    std::vector<std::complex<double> > evens;
    std::vector<std::complex<double> > odds;
    for (int i = 0; i < n; ++i) {
      if (i % 2 == 0) {
        evens.push_back(poly[i]);
      } else {
        odds.push_back(poly[i]);
      }
    }

    // Recurse
    auto solvedE = fftRecursive(evens, power * 2, useOmega, omegas);
    auto solvedO = fftRecursive(odds, power * 2, useOmega, omegas);

    // Construct Results
   std::vector<std::complex<double> > result;
    for (int i = 0; i < n/2; ++i) {
      if (useOmega) {
        result.push_back(solvedE[i] + (omegas.first[i * power] * solvedO[i] ) );
      } else {
        result.push_back(solvedE[i] + (omegas.second[i * power] * solvedO[i] ) );
      }
    }
    for (int i = 0; i < n/2; ++i) {
      if (useOmega) {
        result.push_back(solvedE[i] - (omegas.first[i * power] * solvedO[i] ) );
      } else {
        result.push_back(solvedE[i] - (omegas.second[i * power] * solvedO[i] ) );
      }
    }
    return result;
  }

}

#endif