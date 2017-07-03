#include <stdio.h>
#include <math.h>
#include <vector>
#include <emscripten/emscripten.h>

#define PI 3.14159265

extern "C" {
  int main() {
    printf("Hello world!\n");
    return 0;
  }

/*
 * pre emphasis filter and hamming window
 *
 */
EMSCRIPTEN_KEEPALIVE
void preEmphHamming(std::vector<float>& _signal) {

  // 0.54 - 0.46 * cos( 2.0 * PI * i / (N - 1) ) = 0.08 if N == 1
  // y.at(0) = _signal[0] * 0.08;
  _signal.at(0) *= 0.08;

  for (int i = 1; i < _signal.size(); i++) {
    _signal.at(i) -= 0.97 * _signal.at(i - 1);
    _signal.at(i) *= 0.54 - 0.46 * cos( 2.0 * PI * i / (_signal.size() - 1) );
  }
}

EMSCRIPTEN_KEEPALIVE
void mfcc(float * _signal, size_t length) {
  // TODO: もっといい感じの書き方
  std::vector<float> s(length);
  for(int i=0; i<length; i++) s.at(i) = _signal[i];

  preEmphHamming(s);
  // powerSpectrum(s);
  // lmfb(s);
  // dct(s);
}

}
