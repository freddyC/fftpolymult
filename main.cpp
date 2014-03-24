
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
#include <sstream>

std::pair<double, double> solveIt (int length, bool recursive);
std::vector<double> recursivePolyMult (std::vector<double>, std::vector<double>, int n, std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > &omegas);
// std::vector<double> dynamicPolyMult (std::vector<double>, std::vector<double>, int n, std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > &omegas);
void padPoly (std::vector<double> &poly, std::vector< std::complex <double> > &complex);
std::string writeVector (std::vector <double> );
std::string writeVector (std::vector <std::complex <double> >);
double randomDouble();


int main () {

  std::cout << "\nRecursive Results\n";
  for (long long int i = 2; i < std::pow(2, 20) +1; i *= 2) {
    auto timeTook = solveIt(i, true);
    std::cout << "\n Terms: " << i << " Average: " << timeTook.first << " Deviation: " << timeTook.second << "\n";
  }

  // std::cout << "\nDynamic Results\n";
  // for (long long int i = 2; i < std::pow(2, 20) +1; i *= 2) {
  //   auto timeTook = solveIt(i, false);
  //   std::cout << "\n Terms: " << i << " Average: " << timeTook.first << " Deviation: " << timeTook.second << "\n";
  // }

  return 0;
}

std::pair<double, double> solveIt(int size, bool recursive) {
  std::vector<double> x;
  std::vector<double> y;
  auto omegas = fft::getOmegasForSize(size * 2);

  // if (recursive) {
    return tools::funcEval([&] () {
      for (int i = 0; i < size; ++i) {
        x.push_back(randomDouble());
        y.push_back(randomDouble());
      }
      recursivePolyMult(x, y, size, omegas);
    }, 5);
  // }
  // return tools::funcEval([&] () {
  //   for (int i = 0; i < size; ++i) {
  //     x[i] = randomDouble();
  //     y[i] = randomDouble();
  //   }
  //   dynamicPolyMult(x, y, size)
  // });
}


std::vector<double> recursivePolyMult(std::vector<double> polyA, std::vector<double> polyB, int n, std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > &omegas) {
  /// Double the size and pad with 0's
  std::vector<std::complex<double> > complexA;
  std::vector<std::complex<double> > complexB;
  padPoly(polyA, complexA);
  padPoly(polyB, complexB);

  // fftRecursive each polynomial
  auto solveA = fft::fftRecursive(complexA, 1, true, omegas);
  auto solveB = fft::fftRecursive(complexB, 1, true, omegas);

  // multiply the evaluated polynomials
  std::vector<std::complex<double> > res;
  for (int i = 0; i < solveA.size(); ++i) {
    res.push_back(solveA[i] * solveB[i]);
  }

  // invers fft
  auto complexResult = fft::fftRecursive(res, 1, false, omegas);

  // divide by 2n
  std::vector<double> result;
  for(int i = 0; i < complexResult.size() - 1; ++i ) {
    result.push_back(complexResult[i].real() / (2 * n));
  }

  // TO SEE RESULTS UNCOMENT ME!!
  // std::cout << "\nResults: \n(" << writeVector(polyA) << ")*(" << writeVector(polyB) << ") = (" << writeVector(result) << ")\n\n";
  return result;
}



// std::vector<double> dynamicPolyMult (std::vector<double>, std::vector<double>, int n, std::pair< std::vector< std::complex<double> >, std::vector< std::complex<double> > > &omegas) {


//   return ;
// }


void padPoly(std::vector<double> &poly, std::vector< std::complex <double> > &complex) {
  complex.clear();
  for (int i = 0; i < poly.size(); ++i) {
    complex.push_back(std::complex<double> (poly[i], 0));
  }
  for (int i = 0; i < poly.size(); ++i) {
    complex.push_back(std::complex<double> (0,0));
  }
}

std::string writeVector (std::vector <double> arr) {
  std::stringstream ss;
  ss << arr[0];
  for (int i = 1; i < arr.size(); ++i) {
    ss << " + " << arr[i] << "x^" << i;
  }
  return ss.str();
}

std::string writeVector (std::vector <std::complex <double> > arr) {
  std::stringstream ss;
  ss << arr[0].real();
  for (int i = 1; i < arr.size(); ++i) {
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