#include "ABFM.h"
// NOT THREADSAFE NAITIVELY
void Adams_Bashforth_2_vec(std::vector<double> *Y_next, std::vector<double> &Y_current, std::vector<double> &df_dts, std::vector<double> &previous_df_dts,  double step_size) {
  // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))
  for (unsigned int i = 0; i < (Y_current).size(); i++) {
    (*Y_next)[i] = (Y_current)[i] + (step_size / 2) * ((3 * (df_dts)[i]) - (previous_df_dts)[i]);
  }
  return;
}
void Adams_Bashforth_2(double *Y_next,  double &Y_current,  double &df_dts,  double &previous_df_dts,  double step_size) {
  // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))
    (*Y_next) = (Y_current) + (step_size / 2) * ((3 * df_dts) - (previous_df_dts));
  return;
}

void Forward_Euler_vec(std::vector<double> *Y_next,  std::vector<double> &Y_current,  std::vector<double> &df_dts,  double step_size) {
  // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))
  for (unsigned int i = 0; i < (Y_current).size(); i++) {
    (*Y_next)[i] = (Y_current)[i] + step_size*df_dts[i];
  }
  return;
}
void Forward_Euler(double *Y_next,  double &Y_current,  double &df_dts,  double step_size) {
  // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))
    (*Y_next) = (Y_current) + step_size*df_dts;
  return;
}
