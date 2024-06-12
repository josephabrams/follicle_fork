#include "./conversions.h"
/*convert from molarity to molality using the polynomial Ax^3+Bx^2+Cx+D, D is 0 where they intercept 
 * from python files in virial folder using CRC data located there for 20 degrees C --assumed ~= at 23 C
 * To convert from  molarity  to  molality  for  EG  the coefficients are:  [ 0.01277023 -0.02080002  1.13628309  0.        ]
 * To convert from  molarity  to  molality  for  GLY  the coefficients are:  [0.00665793 0.06926148 1.00285049 0.        ]
 * To convert from  molarity  to  molality  for  NaCl  the coefficients are:  [ 5.80771883e-05 -1.93492494e-02  1.00003795e+00  0.00000000e+00]
 * */
std::unordered_map<std::string,double> molal_conversion_coeff_A={{"EG",0.01277023},{"GLY",0.00665793},{"NaCl",5.80771883e-05}};
std::unordered_map<std::string,double> molal_conversion_coeff_B={{"EG",-0.02080002},{"GLY",0.06926148},{"NaCl",-1.93492494e-02}};
std::unordered_map<std::string,double> molal_conversion_coeff_C={{"EG",1.13628309},{"GLY",1.00285049},{"NaCl",1.00003795}};
std::unordered_map<std::string,double> molal_conversion_coeff_D={{"EG",0.0},{"GLY",0.0},{"NaCl",0.0}};

double molarity_to_molality(double molarity, std::string component_name)
{
  double molality=0.0;
  molality= molal_conversion_coeff_A[component_name]*molarity*molarity*molarity+molal_conversion_coeff_B[component_name]*molarity*molarity+molal_conversion_coeff_C[component_name]*molarity+molal_conversion_coeff_D[component_name];
  return molality;
}
/*
 * The virial osmotic coefficients for the cubic virial osmotic equation and molar_mass
 * */
std::unordered_map<std::string,double> virial_coeff_B={{"EG",0.037},{"GLY",0.023},{"NaCl",0.044}};
std::unordered_map<std::string,double> virial_coeff_C={{"EG",-0.001},{"GLY",0.0},{"NaCl",0.0}};
std::unordered_map<std::string,double> kdiss={{"EG",1.0},{"GLY",1.0},{"NaCl",1.678}};
std::unordered_map<std::string,double> molar_mass={{"EG",62.07},{"GLY",92.09},{"NaCl",58.44}};

double binary_virial(double molality, std::string component_name)
{//osmolality=mi+B_i*m_i^2+C_i*m_i^3
  double osmolality=molality+virial_coeff_B[component_name]*molality*molality+virial_coeff_C[component_name]*molality*molality*molality;
  return osmolality;
}

/* The following equation computes the osmolality of a ternary solution with m ininitial molalilities using B and C coefficients,with kdiss for ions*/


double ternary_virial(double molality_1, double molality_2, std::string component_1, std::string component_2)
{
    //compute virial equation uses vectors for future expansion
  double osmolality=0.0;
  // #pragma omp reduction(+:osmolality) 
  // avoid thread issues by mearly computing all the terms
  std::vector <double> m={molality_1*kdiss[component_1],molality_2*kdiss[component_2]};
  std::vector <double> B={virial_coeff_B[component_1],virial_coeff_B[component_2]};
  std::vector <double> C={virial_coeff_C[component_1],virial_coeff_C[component_2]};
  osmolality=m[0]+m[1]+B[0]*m[0]*m[0]+B[1]*m[1]*m[1]+(B[0]+B[1])*m[0]*m[1]+C[0]*m[0]*m[0]*m[0]+C[1]*m[1]*m[1]*m[1]+3*(std::pow((C[0]*C[0]*C[1]),(1.0/3.0))*m[1]*m[1]*m[2])+3*(std::pow((C[0]*C[1]*C[1]),(1.0/3.0))*m[1]*m[2]*m[2]);
  // for (size_t i = 0; i < m.size(); i++)
  // {
  //   osmolality+=m[i];
  //   for (size_t j = 0; j < B.size(); j++)
  //   {
  //       osmolality+=((B[i]+B[j])/2)*m[i]*m[j];
  //       for (size_t k = 0; k < C.size(); k++)
  //       {
  //         std::cout<<"Cs"<< std::abs(C[i]*C[j]*C[k])<<"\n";
  //         std::cout<<"test pow "<< (std::pow(std::abs(C[i]*C[j]*C[k]),(0.333333)))<<"\n";
  //         std::cout<<"test ms: "<<m[i]*m[j]*m[k]<<"\n";
  //           osmolality+=(std::pow((C[i]*C[j]*C[k]),(0.333333)))*m[i]*m[j]*m[k];
  //       }
  //       
  //   }
  // }
  // std::cout<<"Osmole "<< osmolality<< "\n";
    return osmolality;
}


std::unordered_map<std::string,double> solute_specific_volume={{"NaCl", 0.01661}, {"EG", 0.0557414}, {"GLY", 0.0730903}, {"HM", 0.01661}};//{HM,EG,GLY,PBS} in um^3/femptomole (l/mole)
//eventually pull values from xml for solution thermodynamic constants
// const std::vector <double> Spring_Cell::solute_specific_volume={0.0,0.01661,0.0557414,0.0730903,0.01661};//{0,HM,EG,GLY,PBS} in um^3/femptomole (l/mole)
//
double moles_to_volume(double moles, std::string component_name)
{
  double volume=0.0;
  volume=solute_specific_volume[component_name]*moles;
  return volume;
}
