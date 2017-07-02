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
 * pre emphasis filter
 *   y(n) = x(n) - p x(n - 1)
 */
EMSCRIPTEN_KEEPALIVE
void preEmphasis(float * _signal, size_t length, float p) {
  std::vector<float> y(length);

  printf("signal at 200 is %f\n", _signal[200]);
  printf("signal at 400 is %f\n", _signal[400]);
  printf("signal at 800 is %f\n", _signal[800]);
  printf("signal at 1600 is %f\n", _signal[1600]);

  y.at(0) = _signal[0];

  for(int i = 1; i < length; i++) {
    y.at(i) = _signal[i] - ( p * _signal[i - 1]);
  }

  printf("filtered signal at 200 is %f\n", y[200]);
  printf("filtered signal at 400 is %f\n", y[400]);
  printf("filtered signal at 800 is %f\n", y[800]);
  printf("filtered signal at 1600 is %f\n", y[1600]);
}

/*
 * hamming window
 *   w(n) = 0.54 - 0.46 cos( 2 PI n / (M - 1) ), 0 <= n <= M - 1
 */
// std::array<float, N> hamming(std::array<float, N> _signal) {
//   std::array<float, N> w;
//   for(int i = 0; i < N; i += 1) {
//     w[i] = _signal[i] * (0.54 - 0.46 * cos( 2.0 * PI * i / (N - 1)));
//   }
//
//   return w;
// }

}
