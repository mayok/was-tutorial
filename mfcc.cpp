#include <complex>
#include <vector>
#include <map>
#include <math.h>

class MFCC {
private:
  const double PI = 4 * atan(1.0);
  std::vector<std::complex<double>> Xjo, Xjo2;
  std::map<int, std::map<int, std::complex<double>>> bitTable_;
  size_t numFFT;

private:
  inline double hz2mel(double f) {
    return 2595 * std::log10(1 + f/700);
  }

  inline double mel2hz(double m) {
    return 700 * (std::pow(10, m/2595) - 1);
  }

// constructor から呼び出されるものたち
private:
  // fft で使う
  void bitTable(void) {
    const std::complex<double> j(0, 1); // 0 + j
    for(int n = 2; n <= numFFT; n *= 2) {
      for(int k = 0; k <= N/2 - 1; k+=1) {
        bitTable_[N][k] = exp(-2 * PI * k / n * j);
      }
    }
  }

private:
  // preemphasis filter をかける
  // hamming window をかける
  // fft
  // 振幅スペクトルを求める
  // log mel filter bank をかける
  // 離散コサイン変換 (dct2)
  void preEmphasisAndHamming(void) {
    std::vector<double> f(frame.size());
    f[0] = frame[0] * 0.08;
    for(int i = 1; i < frame.size(); i+=1) {
      // TODO: preEmphasis Coef
      f[i] = frame[i] - ( 0.97 * frame[i - 1]);
      // TODO: N
      f[i] = f[i] * (0.54 - 0.46 * cos( 2.0 * PI * i / (N - 1)));
    }
    frame = f;
  }

  std::vector<std::complex<double>> fft(std::vector<std::complex<double>> x) {
    int N = x.size();
    if(N==1) return x;

    std::vector<std::complex<double>> xe(N/2, 0), xo(N/2, 0);
    for(int i = 0; i < N; i += 2)
      xe[i/2] = x[i];
    for(int i = 1; i < N; i += 2)
      xo[(i-1)/2] = x[i];

    Xjo = fft(xe);
    Xjo2 = fft(xo);
    Xjo.insert(Xjo.end(), Xjo2.begin(), Xjo2.end());

    for(int i = 0; i <= N/2 - 1; i+=1) {
      std::complex<double> t = Xjo[i], tw = bitTable_[N][i];
      Xjo[i] = t + tw * Xjo[i+N/2];
      Xjo[i+N/2] = t - tw * Xjo[i+N/2];
    }

    return Xjo;
  }

  void powerSpectrum(void) {
    // std::vector<std::complex<double>> data
    // for(int i = 0; i < ; i+=1)
    // pow(abs(data[i]),2)
  }

  void logMelFilterbank(void) {
    // いろいろ
  }

  /*
   * DCT: Discrete Cosine Transform
   *
   *            N-1
   *  y[k] = 2* sum x[n]*cos(pi*k*(2n+1)/(2*N)), 0 <= k < N.
   *            n=0
   */
  void dct(void) {

  }

public:
  // constructor
  MFCC() {

    numFFT = 2048;
    bitTable();
  }

  // main process
  // TODO: should be return some value
  void process() {
    // frame = wav.getChannelData(0);
    preEmphasisAndHamming();
    fft();
    // power spectrum ( or amplitude spectrum )
    // apply log mel filter bank
    // dct
  }
}
