#include <stdio.h>
#include <math.h>
#include <array>
#include <emscripten/emscripten.h>

#define N 8192
#define PI 3.14159265

#ifdef __cplusplus
extern "C" {
#endif
/*
 * pre emphasis filter
 *   y(n) = x(n) - p x(n - 1)
 */
EMSCRIPTEN_KEEPALIVE
char * preEmphasis(std::array<float, N> _signal[N], float p) {
  static char y[N];
  y[0] = _signal[0];

  for(int i = 0; i < N; i++) {
    y[i] = _signal[i] - ( p * _signal[i - 1]);
  }

  return y;
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

#ifdef __cplusplus
}
#endif
