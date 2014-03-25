#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex>
#include <cmath>
#include <utility>
#include <iostream>

namespace fft
{
  const double PI = 3.14159265359;
  
  std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > getOmegasForSize (int n) {
    std::vector<std::complex<double> > omega(n);
    std::vector<std::complex<double> > invOmega(n);
    
    for (auto i = 0; i < n; ++i) {
      omega[i] = std::complex<double> (std::cos((2 * PI * i)/n), std::sin((2 * PI * i)/n));
      invOmega[i] = std::complex<double> (std::cos((2 * PI * i)/n), -1 * std::sin((2 * PI * i)/n));
    }
    return std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > (omega, invOmega);
  }
  
  
  int rbs(int i, int k) {
    if (k == 0) {
      return i;
    } else if (i %2 == 1) {
      return std::pow(2, k - 1) + rbs(i / 2, k - 1);
    }
    return rbs(i / 2, k - 1);
  }
  
  std::vector<std::vector<int> > getRbsCache (int size) {
    int logN = std::log(size) / std::log(2);
    std::vector<std::vector<int> > res(size, std::vector<int> (logN +1));
    for (auto i = 0; i < size; ++i) {
      for (auto j = 0; j <= logN; ++j) {
        res[i][j] = rbs(i, logN);
      }
    }
    return res;
  }
  
  std::vector<std::complex<double> > fftRecursive(std::vector<std::complex<double> > &poly, long long int power, std::vector< std::complex<double> > &omegas) {
      // BASE CASE
    if (poly.size() == 1) {
      return poly;
    }
    
      // Divide into two parts
    int n = poly.size();
    std::vector<std::complex<double> > evens;
    std::vector<std::complex<double> > odds;
    for (auto i = 0; i < n; ++i) {
      if (i %2 == 0) {
        evens.push_back(poly[i]);
      } else {
        odds.push_back(poly[i]);
      }
    }
    
      // Recurse
    auto solvedE = fftRecursive(evens, power * 2, omegas);
    auto solvedO = fftRecursive(odds, power * 2, omegas);
    
      // Construct Results
    std::vector<std::complex<double> > result;
    std::vector<std::complex<double> > right;
    
    for (auto i = 0; i < n/2; ++i) {
      right.push_back(omegas[i * power] * solvedO[i] );
    }
    
    for (auto i = 0; i < n/2; ++i) {
      result.push_back(solvedE[i] + right[i]);
    }
    for (auto i = 0; i < n/2; ++i) {
      result.push_back(solvedE[i] - right[i]);
    }
    return result;
  }
  
  std::vector<std::complex<double> > fftDynamic (std::vector<std::complex<double> > &poly, std::vector< std::complex<double> > &omegas, std::vector<std::vector<int> > &rbsCache) {
    int n = poly.size();
    int logN = std::log(n) / std::log(2);
    
      // set up 2d vector of solutions
    std::vector<std::vector<std::complex<double> > > sol(2, std::vector<std::complex<double> > (n, std::complex<double> (0,0)));
    
      // Add in bottom row (base case)
    for (auto i = 0; i < n; ++i) {
      sol[0][rbsCache[i][logN]] = poly[i];
    }
    
      // number of omega slots at the bottom
    int power = n/2;
    int size = 2;
    
      // scan from bottom to top
    for (auto k = 1; k <= logN; ++k) {
        // scan accross by size for each sub solution
      for (auto i = 0; i < n; i += size) {
        std::complex<double> odd;
          // fill in this solution
        for (auto j = 0; j < size/2; ++j) {
          odd = omegas[j*power] * sol[(k-1) %2][i+j+(size/2)];
          sol[k %2][i+j] = sol[(k-1) %2][i+j] + odd;
          sol[k %2][i+j+(size/2)] = sol[(k-1) %2][i+j] - odd;
        }
      }
        // decrease power and increase size
      power /= 2;
      size *=2;
    }
    return sol[logN %2];
  }
}

#endif