#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#define WITHOUT_NUMPY
#include "matplotlibcpp.h"

#define PI 4 * atan(1.0)

namespace plt = matplotlibcpp;

void preEmphasisAndHamming(std::vector<double> frame) {
  int N = frame.size();
  static std::vector<double> f(N);

  frame[0] *= 0.08;
  for(int i = 1; i < N; i++) {
    f[i] = frame[i] - ( 0.97 * frame[i - 1]);
    frame[i] = f[i] * (0.54 - 0.46 * cos( 2.0 * PI * i / (N - 1)));
  }

  // pre-emphasis and hamming filtered signal
  plt::plot(frame);
  plt::show();
}

void process(int16_t* samples, size_t N) {
  std::vector<double> frame(N);

  for(int i = 0; i < N; i++) {
    frame.at(i) = samples[i];
  }

  // original signal
  plt::plot(frame);
  plt::show();

  preEmphasisAndHamming(frame);
}

int main(int argc, char ** argv) {
  if(argc < 1) return 1;
  std::ifstream wavFp;
  wavFp.open(argv[1]);
  if(!wavFp.is_open()) return 1;

  uint16_t bufferLength = 1764;
  int16_t* buffer = new int16_t[bufferLength];
  int bufferBPS = (sizeof buffer[0]);

  wavFp.seekg(6026);
  wavFp.read((char *)buffer, bufferLength*bufferBPS);

  process(buffer, bufferLength);

  wavFp.close();
}
