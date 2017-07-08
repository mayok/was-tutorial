#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>

using namespace Eigen;

// ガウス分布の数
#define K 2

// 次元数
#define D 2

// データの数
#define N 100

#define PI 4*atan(1.0)

typedef Eigen::Matrix<float, 1, D> m_d;

// 多次元 (多変量) ガウス分布
float gaussian(VectorXd x, VectorXd mu, MatrixXd sigma) {
  return exp( -0.5 * (x - mu).dot( sigma.inverse() * (x - mu).transpose() ) )
    / pow(sqrt(2 * PI), D) * sqrt(sigma.determinant());
}

int main() {
  // N 行 D 列
  std::vector<Matirx<float, N, D> x(N, Matrix<float, N, D>::Random());

  // 平均 mu, 分散 sigma, 混合係数 pi を初期化する
  // 平均
  /*
      [
        [mu_x, mu_y], // 次元数 (D) 個
        [mu_x, mu_y],
        ... K 個
      ]
      ガウス分布 K に対する, X軸の平均，Y軸の平均 ...
   */
  std::vector<m_d> mu(K, m_d::Random());

  // 分散 (分散共分散行列)
  // D x D 行列
  std::vector<Matrix<float,D,D> > sigma(K, Matrix<float,D,D>::Random());

  // TODO: 制約 Sigma(k) pi_k = 1 を満たすように初期化する
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  std::vector<float> pi(K);
  for(int i=0; i<K; i++) {
    pi[i] = distribution(mt);
  }

  // TODO: gamma, N, K

  while(true) {
    // E-step: パラメータ (mu, sigma, pi) を使って負担率 gamma を計算する
    for(int n = 0; n < N; n++) {
      float t = 0.0;
      for(int k = 0; k < K; k++) {
        t += pi[k] * gaussian(x[n], mu[k], sigma[k]);
      }

      for(int k = 0; k < K; k++)
        gamma[n][k] = pi[k] * gaussian(x[n], mu[k], sigma[k]) / t;
    }

    // M-step: 負担率を使ってパラメータを更新する
    for(int k = 0; k < K; k++) {
      // N_k
      float Nk = 0.0;
      for(int n = 0; n < N; n++) {
        Nk += gamma[n][k];
      }

      // 平均
      float _mu = 0.0;
      for(int n = 0; n < N; n++) {
        _mu += gamma[n][k] * x[n]
      }
      mu[k] = _mu / Nk;

      // 分散
      MatrixXd _sigma;
      for(int n = 0; n < N; n++) {
        _sigma += gamma[n][k] * (x[n] - mu[k]) * (x[n] - mu[k]).transpose();
      }
      sigma[k] = _sigma / Nk;

      // 混合係数
      pi[k] = Nk / N;

    }

    // 収束性の確認
    // 対数尤度関数 Sigma(n) { log ( Sigma(k) pi * N )}

    // if 収束 break;
  }
}
