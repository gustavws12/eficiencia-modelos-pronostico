// =========================================================================
// hamilton.cpp — Filtro de Hamilton en C++ para el MS-AR
//
// La verosimilitud del MS-AR se evalúa miles de veces dentro del optimizador
// y decenas de miles de veces en el Monte Carlo. En R puro es el cuello de
// botella del experimento (≈ 8 s por ajuste); aquí baja a milisegundos.
// Si Rcpp no está disponible, 02_modelos.R usa la versión en R puro.
// =========================================================================
#include <Rcpp.h>
using namespace Rcpp;

// Distribución ergódica de P por iteración de potencia (K pequeño)
static NumericVector ergodica(const NumericMatrix& P, int K) {
  NumericVector xi(K, 1.0 / K), nx(K);
  for (int it = 0; it < 300; ++it) {
    for (int j = 0; j < K; ++j) {
      double acc = 0.0;
      for (int i = 0; i < K; ++i) acc += P(i, j) * xi[i];
      nx[j] = acc;
    }
    double s = 0.0;
    for (int j = 0; j < K; ++j) s += nx[j];
    if (!(s > 0)) { for (int j = 0; j < K; ++j) xi[j] = 1.0 / K; return xi; }
    for (int j = 0; j < K; ++j) xi[j] = nx[j] / s;
  }
  return xi;
}

// Verosimilitud negativa directamente desde el vector de parámetros.
// Evita todo el trabajo de R dentro del optimizador (desempaque, softmax,
// distribución ergódica), que era el costo dominante del MS-AR.
// [[Rcpp::export]]
double msar_nll_cpp(NumericVector theta, NumericVector Y, NumericMatrix X,
                    int K, int orden) {
  int n = Y.size();
  NumericVector cs(K), sig(K);
  NumericMatrix phis(K, orden), P(K, K);

  for (int k = 0; k < K; ++k) cs[k] = theta[k];
  for (int j = 0; j < orden; ++j)
    for (int k = 0; k < K; ++k) phis(k, j) = theta[K + j * K + k];
  for (int k = 0; k < K; ++k) {
    double ls = theta[K + K * orden + k];
    if (ls < -6.0) ls = -6.0;
    if (ls >  3.0) ls =  3.0;
    sig[k] = std::exp(ls);
  }
  int pos = 2 * K + K * orden;
  for (int i = 0; i < K; ++i) {
    NumericVector fila(K);
    for (int j = 0; j < K; ++j) fila[j] = (j == i) ? 0.0 : theta[pos++];
    double s = 0.0;
    for (int j = 0; j < K; ++j) {
      double v = fila[j];
      if (v < -20.0) v = -20.0;
      if (v >  20.0) v =  20.0;
      fila[j] = std::exp(v); s += fila[j];
    }
    for (int j = 0; j < K; ++j) P(i, j) = fila[j] / s;
  }

  NumericVector xi = ergodica(P, K);
  const double LOG_SQRT_2PI = 0.9189385332046727;
  double ll = 0.0;
  NumericVector pred(K), num(K);

  for (int t = 0; t < n; ++t) {
    for (int j = 0; j < K; ++j) {
      double acc = 0.0;
      for (int i = 0; i < K; ++i) acc += P(i, j) * xi[i];
      pred[j] = acc;
    }
    double f = 0.0;
    for (int j = 0; j < K; ++j) {
      double mu = cs[j];
      for (int q = 0; q < orden; ++q) mu += phis(j, q) * X(t, q);
      double zs = (Y[t] - mu) / sig[j];
      double dens = std::exp(-0.5 * zs * zs - LOG_SQRT_2PI) / sig[j];
      if (!R_finite(dens) || dens < 1e-300) dens = 1e-300;
      num[j] = pred[j] * dens;
      f += num[j];
    }
    if (!R_finite(f) || f <= 0.0) return 1e10;
    ll += std::log(f);
    for (int j = 0; j < K; ++j) xi[j] = num[j] / f;
  }
  if (!R_finite(ll)) return 1e10;
  return -ll;
}

// [[Rcpp::export]]
List hamilton_filtro(NumericVector Y, NumericMatrix X,
                     NumericVector cs, NumericMatrix phis,
                     NumericVector sig, NumericMatrix P,
                     NumericVector xi0, bool guardar_filtro) {
  int n = Y.size();
  int K = cs.size();
  int p = X.ncol();

  NumericMatrix eta(n, K);
  const double LOG_SQRT_2PI = 0.9189385332046727;
  for (int k = 0; k < K; ++k) {
    double s = sig[k];
    if (!(s > 0) || !R_finite(s)) return List::create(_["ok"] = false);
    for (int t = 0; t < n; ++t) {
      double mu = cs[k];
      for (int j = 0; j < p; ++j) mu += phis(k, j) * X(t, j);
      double zscore = (Y[t] - mu) / s;
      double dens = std::exp(-0.5 * zscore * zscore - LOG_SQRT_2PI) / s;
      if (!R_finite(dens) || dens < 1e-300) dens = 1e-300;
      eta(t, k) = dens;
    }
  }

  NumericVector xi = clone(xi0);
  NumericVector pred(K), num(K);
  NumericMatrix filt(guardar_filtro ? n : 1, K);
  double ll = 0.0;

  for (int t = 0; t < n; ++t) {
    for (int j = 0; j < K; ++j) {
      double acc = 0.0;
      for (int i = 0; i < K; ++i) acc += P(i, j) * xi[i];
      pred[j] = acc;
    }
    double f = 0.0;
    for (int j = 0; j < K; ++j) { num[j] = pred[j] * eta(t, j); f += num[j]; }
    if (!R_finite(f) || f <= 0.0) return List::create(_["ok"] = false);
    ll += std::log(f);
    for (int j = 0; j < K; ++j) xi[j] = num[j] / f;
    if (guardar_filtro) for (int j = 0; j < K; ++j) filt(t, j) = xi[j];
  }

  return List::create(_["ok"] = true, _["ll"] = ll, _["xi_T"] = xi,
                      _["filt"] = filt);
}
