
#include "tools.h"
#include "fft.h"
#include <vector>
#include <utility>
#include <iostream>
#include <random>
#include <algorithm>
#include <complex>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

std::pair<double, double> solveIt (long long int length, bool recursive);
std::vector<double> recursivePolyMult (std::vector<double> &polyA, std::vector<double> &polyB, int n, std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > &omegas);
std::vector<double> dynamicPolyMult (std::vector<double> &polyA, std::vector<double> &polyB, int n, std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > &omegas, std::vector<std::vector<int> > &rbsCache);
void padPoly (std::vector<double> &poly, std::vector< std::complex <double> > &complex);
std::string writeVector (std::vector <double> );
std::string writeVector (std::vector <std::complex <double> >);
double randomDouble();


int main () {
    // Do the recursive
    // std::ofstream foutR ("/Users/fred.christensen/Dropbox/school/Algorithms/fftpolymult/recursive.csv");
    // foutR << "Polynomial Length, Time Took" << std::endl;
    // std::cout << "\nRecursive Results\n";
    // for (auto i = 2; i < std::pow(2, 20) +1; i *= 2) {
    //   auto timeTook = solveIt(i, true);
    //   std::cout << "\n Terms: " << i << " Average: " << timeTook.first << "\n";
    // foutR << i << ", " << timeTook.first << std::endl;
    // }
    // foutR.close();
  
    // do the dynamic
  std::ofstream foutD ("/Users/fred.christensen/Dropbox/school/Algorithms/fftpolymult/dynamic.csv");
  foutD << "Polynomial Length, Time Took" << std::endl;
  std::cout << "\nDynamic Results\n";
  for (auto i = 2; i < std::pow(2, 20) +1; i *= 2) {
    auto timeTook = solveIt(i, false);
    foutD << i << ", " << timeTook.first << std::endl;
    std::cout << "\n Terms: " << i << " Average: " << timeTook.first << "\n";
  }
  foutD.close();
  return 0;
}

std::pair<double, double> solveIt(long long int size, bool recursive) {
  std::vector<double> x;
  std::vector<double> y;
  auto omegas = fft::getOmegasForSize(size * 2);
  
  if (recursive) {
    return tools::funcEval([&] () {
      x.clear();
      y.clear();
      for (auto i = 0; i < size; ++i) {
        x.push_back(randomDouble());
        y.push_back(randomDouble());
      }
      recursivePolyMult(x, y, size, omegas);
    });
  }
  
  auto rbsCache = fft::getRbsCache(size * 2);
  return tools::funcEval([&] () {
    x.clear();
    y.clear();
    for (auto i = 0; i < size; ++i) {
      x.push_back(randomDouble());
      y.push_back(randomDouble());
    }
    dynamicPolyMult(x, y, size, omegas, rbsCache);
  });
}


std::vector<double> recursivePolyMult (std::vector<double> &polyA, std::vector<double> &polyB, int n, std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > &omegas) {
    /// Double the size and pad with 0's
  std::vector<std::complex<double> > complexA;
  std::vector<std::complex<double> > complexB;
  padPoly(polyA, complexA);
  padPoly(polyB, complexB);
  
    // fftRecursive each polynomial
  auto solveA = fft::fftRecursive(complexA, 1, omegas.first);
  auto solveB = fft::fftRecursive(complexB, 1, omegas.first);
  
    // multiply the evaluated polynomials
  std::vector<std::complex<double> > res;
  for (auto i = 0; i < solveA.size(); ++i) {
    res.push_back(solveA[i] * solveB[i]);
  }
  
    // invers fft
  auto complexResult = fft::fftRecursive(res, 1, omegas.second);
  
    // divide by 2n
  std::vector<double> result;
  for (auto i = 0; i < complexResult.size() - 1; ++i ) {
    result.push_back(complexResult[i].real() / (2 * n));
  }
  
    // TO SEE RESULTS UNCOMENT ME!!
    // std::cout << "\nResults: \n(" << writeVector(polyA) << ")*(" << writeVector(polyB) << ") = (" << writeVector(result) << ")\n\n";
  return result;
}



std::vector<double> dynamicPolyMult (std::vector<double> &polyA, std::vector<double> &polyB, int n, std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > &omegas, std::vector<std::vector<int> > &rbsCache) {
    // Double the size and pad with 0's
  std::vector<std::complex<double> > complexA;
  std::vector<std::complex<double> > complexB;
  padPoly(polyA, complexA);
  padPoly(polyB, complexB);
  
    // fftRecursive each polynomial
  auto solveA = fft::fftDynamic(complexA, omegas.first, rbsCache);
  auto solveB = fft::fftDynamic(complexB, omegas.first, rbsCache);
  
    // multiply the evaluated polynomials
  std::vector<std::complex<double> > res;
  for (auto i = 0; i < solveA.size(); ++i) {
    res.push_back(solveA[i] * solveB[i]);
  }
  
    // invers fft
  auto complexResult = fft::fftDynamic(res, omegas.second, rbsCache);
  
    // divide by 2n
  std::vector<double> result;
  for (auto i = 0; i < complexResult.size() - 1; ++i ) {
    result.push_back(complexResult[i].real() / (2 * n));
  }
  
    // TO SEE RESULTS UNCOMENT ME!!
    // std::cout << "\nResults: \n(" << writeVector(polyA) << ")*(" << writeVector(polyB) << ") = (" << writeVector(result) << ")\n\n";
  return result;
}


void padPoly(std::vector<double> &poly, std::vector< std::complex <double> > &complex) {
  complex.clear();
  for (auto i = 0; i < poly.size(); ++i) {
    complex.push_back(std::complex<double> (poly[i], 0));
  }
  for (auto i = 0; i < poly.size(); ++i) {
    complex.push_back(std::complex<double> (0,0));
  }
}

std::string writeVector (std::vector <double> arr) {
  std::stringstream ss;
  ss << arr[0];
  for (auto i = 1; i < arr.size(); ++i) {
    ss << " + " << arr[i] << "x^" << i;
  }
  return ss.str();
}

std::string writeVector (std::vector <std::complex <double> > arr) {
  std::stringstream ss;
  ss << arr[0].real();
  for (auto i = 1; i < arr.size(); ++i) {
    ss << " + " << arr[i].real() << "x^" << i;
  }
  return ss.str();
}



double randomDouble() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-10, 10);
  return dis(gen);
}