#include <vector>
#include <cmath>

// ガウス分布の数
#define K 2

// 次元数
#define D 2

// データの数
#define N 100

float gaussian(x, mu, sigma) {
  // 多次元 (多変量) ガウス分布
}

int main() {
  // 平均 mu, 分散 sigma, 混合係数 pi を初期化する
  std::vector<float> mu(K);

  // TODO: sigma は行列
  std::vector<float> sigma(K);

  // TODO: 制約 Sigma(k) pi_k = 1 を満たすように初期化する
  std::vector<float> pi(K);

  while(true) {
    // E-step: パラメータ (mu, sigma, pi) を使って負担率 gamma を計算する
    for(int n = 0; n < N; n++) {
      float t = 0.0;
      for(int k = 0; k < K; k++) {
        t += pi[k] * gaussian(x[n], mu[k], sigma[k]);
      }

      for(int k = 0; k < K; k++)
        gaussian[n][k] = pi[k] * gaussian(x[n], mu[k], sigma[k]) / t;
    }

    // M-step: 負担率を使ってパラメータを更新する
    for(int k = 0; k < K; k++) {
      // N_k
      float Nk = 0.0;
      for(int n = 0; n < N; n++) {
        Nk += gamma[n][k];
      }

      // 平均
      float t = 0.0;
      for(int n = 0; n < N; n++) {
        t = gamma[n][k] * x[n]
      }
      mu[k] = t / Nk;

      // 分散
      float t = 0.0;
      for(int n = 0; n < N; n++) {
        t = gamma[n][k] * (x[n] - mu[k]) * ( (x[n] - mu[k]) の転地行列 )
      }
      sigma[k] = t / Nk;

      // 混合係数
      pi[k] = Nk / N;

    }

    // 収束性の確認
    // 対数尤度関数 Sigma(n) { log ( Sigma(k) pi * N )}

    // if 収束 break;
  }
}
