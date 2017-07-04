#include <stdio.h>
#include <math.h>
#include <vector>
#include <complex>
#include <map>
#include <algorithm>
#include <emscripten/emscripten.h>

#define PI 3.14159265
#define NUM_FFT 2048

extern "C" {
  int main() {
    printf("Hello world!\n");
    return 0;
  }

  EMSCRIPTEN_KEEPALIVE
  std::vector<std::complex<double> > fft(std::vector<std::complex<double> > x) {
    int N = x.size();
    if (N==1) return x;

    std::vector<std::complex<double> > xe(N/2,0), xo(N/2, 0), Xjo, Xjo2;
    int i;

    for(i=0;i<N;i+=2) xe[i/2] = x[i];
    for(i=1;i<N;i+=2) xo[(i-1)/2] = x[i];

    Xjo = fft(xe);
    Xjo2 = fft(xo);
    Xjo.insert (Xjo.end(), Xjo2.begin(), Xjo2.end());

    // TODO: この部分を先に計算しておきたい
    std::map<int,std::map<int,std::complex<double> > > bt;
    const std::complex<double> J(0,1);
    for(int n=2;n<=NUM_FFT;n*=2) {
      for(int k=0;k<=n/2-1;k++) {
        bt[n][k] = exp(-2*PI*k/n*J);
      }
    }

    for(i=0;i<N/2-1;i++) {
      std::complex<double> t = Xjo[i], tw = bt[N][i];
      Xjo[i] = t + tw * Xjo[i+N/2];
      Xjo[i+N/2] = t = tw * Xjo[i+N/2];
    }

    return Xjo;
  }

 /*
  * mel filter bank
  *
  */
  EMSCRIPTEN_KEEPALIVE
  void melFilterBank(float filterbank[20][20]) {
    // int fmax = 44100 / 2; // sampleRate / 2
    float melmax = 1127.01048 * log( 22050 / 700.0 + 1.0);
    int nmax = NUM_FFT / 2;

    float df = 44100 / NUM_FFT;
    float dmel = melmax / 21; // melmax / ( number of channels + 1)

    float melcenters[20], fcenters[20], indexcenters[20];
    for(int i = 0; i < 20; i++) {
      melcenters[i] = (i+1) * dmel;
      fcenters[i] = 700.0 * ( exp(melcenters[i] / 1127.01048) - 1);
      indexcenters[i] = round(fcenters[i] / df);
    }

    float indexstart[20], indexstop[20];
    indexstart[0] = 0;
    std::copy( indexcenters, indexcenters + 19, indexstart + 1);
    std::copy( indexcenters + 1, indexcenters + 20, indexstop);
    indexstop[19] = nmax;
    // TODO: ここまでは先に計算できそう

    for(int i = 0; i < 20; i++) {

      float increment = 1.0 / ( indexcenters[i] - indexstart[i] );
      for(int j = indexstart[i]; j < indexcenters[i]; j++)
        filterbank[i][j] = ( j - indexstart[i] ) * increment;

      float decrement = 1.0 / ( indexstop[i] - indexcenters[i]);
      for(int j = indexcenters[i]; j < indexstop[i]; j++)
        filterbank[i][j] = 1.0 - (j - indexcenters[i]) * decrement;
    }

    // fcenters は plot する時に必要になる plot(x: fcenters, y: filterbank);
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

  /*
   *  fft, power spectrum
   *
   */
  EMSCRIPTEN_KEEPALIVE
  void powerSpectrum(std::vector<float>& _signal) {
    _signal.resize(NUM_FFT);
    std::vector<std::complex<double> > fc(_signal.begin(), _signal.end());
    std::vector<std::complex<double> > fftc = fft(fc);

   // power spectrum
   for(int i = 0; i<NUM_FFT/2+1; i++) {
    // _signal[i] = pow(abs(fftc[i]), 2);
    _signal[i] = abs(fftc[i]);
   }
  }

  EMSCRIPTEN_KEEPALIVE
  void lmfb(std::vector<float>& _signal, float mspec[20]) {
    float filterbank[20][20];
    melFilterBank(filterbank);

    //  std::vector<float> mspec(20);
    for(int i = 0; i < 20; i++) {
      float t = 0.0;
      for(int j = 0; j < 20; j++) {
        t += _signal[i] * filterbank[i][j];
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
  void dct(float mspec[20]) {
    float y[20];
    for(int k = 0; k < 20; k++) {
      if(k == 0)
        float t = sqrt(1 / 80);
      float t = sqrt(1 / 40);

      float s = 0.0;
      for(int n = 0; n < 20; n++) {
        s += mspec[n] * cos( PI * k * ( 2 * n + 1 ) / 40.0);
      }

      y[k] = t * 2 * s;
    }
    mspec = y;
  }


  EMSCRIPTEN_KEEPALIVE
  void mfcc(float * _signal, size_t length) {
    float mspec[20];

    // TODO: 配列から vector に代入する, もっといい感じの書き方
    std::vector<float> s(length);
    for(int i=0; i<length; i++) s.at(i) = _signal[i];

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
