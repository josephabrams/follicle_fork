#ifndef __ABFM_H__
#define __ABFM_H__
#include <iostream>
#include <memory>
#include <stdexcept>
#include "../../core/PhysiCell.h"
#include "../../modules/PhysiCell_standard_modules.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <unordered_map>
#include <vector>
// using namespace BioFVM;
// using namespace PhysiCell;

// NOT THREADSAFE NAITIVELY
void Adams_Bashforth_2_vec(std::vector<double> *Y_next,  std::vector<double> &Y_current,  std::vector<double> &df_dts,  std::vector<double> &previous_df_dts,  double step_size);

// void Adams_Bashforth_2_vec(std::vector<double> *Y_next,  std::vector<double> &Y_current,  std::vector<double> &df_dts,  std::vector<double> &previous_df_dts,  double step_size);
void Adams_Bashforth_2(double *Y_next,  double &Y_current,  double &df_dts,  double &previous_df_dts,  double step_size);


void Forward_Euler_vec(std::vector<double> *Y_next,  std::vector<double> &Y_current,  std::vector<double> &df_dts,  double step_size);

void Forward_Euler(double *Y_next,  double &Y_current,  double &df_dts,  double step_size);

#endif

