#include <cmath>
#include <iostream>
#include <string>
#include <array>
using namespace std;

#define _USE_MATH_DEFINES
#define R_ONE 1
#define R_TWO 2
#define K_TWO 5
#define SCREEN_WIDTH 160
#define SCREEN_HEIGHT 25
#define SCREEN_DIM SCREEN_WIDTH*SCREEN_HEIGHT

const double theta_resolution = 0.07;
const double phi_resolution = 0.02;
const double R1 = R_ONE;
const double R2 = R_TWO;
const double K2 = K_TWO;

const double K1 = SCREEN_WIDTH * K_TWO * 3 / (8 * (R_ONE + R_TWO));

template<typename T>
void PrintArray(T container) {
  for (auto _ : container) {
    cout << _;
  }
  cout << endl;
}

template<typename T>
void PrintMatrix(T container) {
  for (auto _ : container) {
    PrintArray(_);
  }
}

template<typename T>
void FillMatrix(T & matrix, char c) {
  for (size_t _ = 0; _ < matrix.size(); _++) {
    matrix[_].fill(c);
  }
}

void RenderFrame(double A, double B) {
  double cosA = cos(A);
  double cosB = cos(B);
  double sinA = sin(A);
  double sinB = sin(B);

  // normally array are y,x
  // init screen
  array<array<char, SCREEN_WIDTH>, SCREEN_HEIGHT> output;
  array<array<char, SCREEN_WIDTH>, SCREEN_HEIGHT> zbuffer;
  FillMatrix(output, ' ');
  FillMatrix(zbuffer, 0);

  for(double theta = 0.0; theta < 2.0*M_PI; theta += theta_resolution) {
    double cos
  }

}

int main() {
  RenderFrame(0.5f, 0.6f);
  cout << "\u03c0 = " << M_PI << endl;
  return 0;
}
