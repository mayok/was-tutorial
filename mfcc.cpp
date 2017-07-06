#include <stdio.h>
#include <cmath>
#include <vector>
#include <complex>
#include <map>
#include <algorithm>
#include <emscripten/emscripten.h>

#define PI 4*atan(1.0)
#define NUM_FFT 2048

extern "C" {
  int main() {
    printf("Hello world!\n");
    return 0;
  }

  // C-linkage ... というエラーが出るので，書き直した
  // EMSCRIPTEN_KEEPALIVE
  // std::vector<std::complex<double> > fft(std::vector<std::complex<double> > x) {
  //   int N = x.size();
  //   if (N==1) return x;
  //
  //   std::vector<std::complex<double> > xe(N/2,0), xo(N/2, 0), Xjo, Xjo2;
  //   int i;
  //
  //   for(i=0;i<N;i+=2) xe[i/2] = x[i];
  //   for(i=1;i<N;i+=2) xo[(i-1)/2] = x[i];
  //
  //   Xjo = fft(xe);
  //   Xjo2 = fft(xo);
  //   Xjo.insert (Xjo.end(), Xjo2.begin(), Xjo2.end());
  //
  //   // TODO: この部分を先に計算しておきたい
  //   std::map<int,std::map<int,std::complex<double> > > bt;
  //   const std::complex<double> J(0,1);
  //   for(int n=2;n<=NUM_FFT;n*=2) {
  //     for(int k=0;k<=n/2-1;k++) {
  //       bt[n][k] = exp(-2*PI*k/n*J);
  //     }
  //   }
  //
  //   for(i=0;i<N/2-1;i++) {
  //     std::complex<double> t = Xjo[i], tw = bt[N][i];
  //     Xjo[i] = t + tw * Xjo[i+N/2];
  //     Xjo[i+N/2] = t = tw * Xjo[i+N/2];
  //   }
  //
  //   return Xjo;
  // }

  // fft
  EMSCRIPTEN_KEEPALIVE
  void fft(std::vector<float>& _signal, std::vector<float>& im) {
    // std::vector<std::complex<float> > でできる
    size_t n = _signal.size();

    for(int i = 0; i < n; i++) {
      int t = 0;
      for(int j = 0, h = i, k = n; ; h >>= 1) {
        k>>=1;
        if(k == 0) {
          t = j;
          break;
        }

        j = (j << 1) | (h & 1);
      }

      if(t > i) {
        std::swap(_signal[i], _signal[t]);
      }
    }
    for(int hn = 1; hn * 2 <= n; hn *= 2) {
      for(int i = 0; i < n; i+= hn * 2) {
        for(int j = i; j < i + hn; j++) {
          float _cos = cos(PI * (j - i) / hn);
          float _sin = sin(PI * (j - i) / hn);
          float tre = _signal[j+hn] * _cos + im[j+hn] * _sin;
          float tim = -1 * _signal[j+hn] * _sin + im[j+hn] * _cos;

          _signal[j+hn] = _signal[j] - tre;
          im[j+hn] = im[j] - tim;
          _signal[j] += tre;
          im[j] += tim;
        }
      }
    }
  }

 /*
  * mel filter bank
  *
  */
  EMSCRIPTEN_KEEPALIVE
  void melFilterBank(float filterbank[20][1024]) {
    int indexcenters[20] = { 6, 13, 21, 31, 42, 55, 71, 90, 112, 138, 169, 205, 248, 299, 358, 429, 512, 610, 726, 863 };
    int indexstart[20] = { 0, 6, 13, 21, 31, 42, 55, 71, 90, 112, 138, 169, 205, 248, 299, 358, 429, 512, 610, 726 };
    int indexstop[20] = { 13, 21, 31, 42, 55, 71, 90, 112, 138, 169, 205, 248, 299, 358, 429, 512, 610, 726, 863, 1024 };

    for(int i = 0; i < 20; i++) {

      float increment = 1.0 / ( indexcenters[i] - indexstart[i] );
      for(int j = indexstart[i]; j < indexcenters[i]; j++)
        filterbank[i][j] = ( j - indexstart[i] ) * increment;

      float decrement = 1.0 / ( indexstop[i] - indexcenters[i]);
      for(int j = indexcenters[i]; j < indexstop[i]; j++)
        filterbank[i][j] = 1.0 - (j - indexcenters[i]) * decrement;
    }
  }

  /*
  * pre emphasis filter and hamming window
  *
  */
  EMSCRIPTEN_KEEPALIVE
  void preEmphHamming(std::vector<float>& _signal) {
    size_t len = _signal.size();
    std::vector<float> y(len, _signal[0] * 0.08);

    for (int i = 1; i < len; i++) {
      y[i] = _signal[i] - ( 0.97 * _signal[i-1]);
      y[i] *= 0.54 - 0.46 * cos( 2.0 * PI * i / (len - 1) );
    }
    _signal = y;
  }

  /*
   *  fft, power spectrum
   *
   */
  EMSCRIPTEN_KEEPALIVE
  void powerSpectrum(std::vector<float>& _signal) {
    _signal.resize(NUM_FFT, 0);
    std::vector<float> im(_signal.size(), 0);
    fft(_signal, im);

   // power spectrum
   for(int i = 0; i<NUM_FFT/2+1; i++) {
    _signal[i] = _signal[i] * _signal[i] + im[i] * im[i];

    // XXX: we use amplitude spectrum
    _signal[i] = sqrt(_signal[i]);
   }
  }

  EMSCRIPTEN_KEEPALIVE
  void lmfb(std::vector<float>& _signal, std::vector<float>& mspec) {
    float filterbank[20][1024];
    melFilterBank(filterbank);

    for(int i = 0; i < 20; i++) {
      float t = 0.0;
      for(int j = 0; j < 1024; j++) {
        t += _signal[j] * filterbank[i][j];
      }
      mspec[i] = log10( t );
    }
  }

  /*
   * DCT: Discrete Cosine Transform
   *
   *            N-1
   *  y[k] = 2* sum x[n]*cos(pi*k*(2n+1)/(2*N)), 0 <= k < N.
   *            n=0
   */
  EMSCRIPTEN_KEEPALIVE
  void dct(std::vector<float>& mspec) {
    std::vector<float> y(20);
    for(int k = 0; k < 20; k++) {
      float s = 0.0;
      for(int n = 0; n < 20; n++) {
        s += mspec[n] * cos( PI * k * ( 2 * n + 1 ) / 40.0);
      }
      if(k == 0)
        y[k] = sqrt(0.0125) * 2 * s;
      else
        y[k] = sqrt(0.025) * 2 * s;
    }
    mspec = y;
  }


  EMSCRIPTEN_KEEPALIVE
  void mfcc(float * _signal, size_t length) {
    std::vector<float> mspec(20);
    std::vector<float> s(std::begin(_signal), std::end(_signal));

    preEmphHamming(s);
    powerSpectrum(s);
    lmfb(s, mspec);
    dct(mspec);

    printf("mfcc: \n");
    for(int i = 0; i < 12; i++) {
      printf("%d: %f\n", i, mspec[i]);
    }
  }

}
