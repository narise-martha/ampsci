#include "Angular/Wigner369j.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/DiagramRPA0_jL.hpp"
#include "HF/HartreeFock.hpp"
#include "Kion_functions.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/UnitConv_conversions.hpp"
#include "StandardHaloModel.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/ostream.hpp"
#include "qip/Methods.hpp"
#include "qip/omp.hpp"
#include <iostream>
#include <memory>
#include <complex>

// RADIAL INTEGRAL FUNCTION
// int func1 * jL * func2 dr
double radial_int(const std::vector<double> &func1,
                  const std::vector<double> &func2,
                  const std::vector<double> &jL, Grid grid) {
  return NumCalc::integrate(grid.du(), 0, 0, grid.drdu(), func1, func2, jL);
}

// ANGULAR INTAGRAL FUNCTIONS

// Clecbsch-Gordon coefficient but in a different order
// because I couldn't deal with the other notation
double cg(int twola, int twolb, int twoL, int twoma, int twomb, int twoM){
  return Angular::cg_2(twola, twoma, twolb, twomb, twoL, twoM);
}

// Matrix elements of a spherical harmonic, Y_{L,M}
// < lb, mb | Y_{L,M} | la, ma >
double ME_spherical_harm(int L, int twoM, int lb, int la, int twomb, int twoma) {
    
    // Reduced Matrix element
    // < lb || C^L_M || la >
    reduced_ME = pow(-1,lb)*sqrt((2*lb+1)*(2*la+1))*Angular::threej_2(2*lb,2*L,2*la,0,0,0);

    // Y_{L,M} = sqrt((2L+1)/4pi) C^L_M
    return sqrt((2*L+1)/4*M_PI)*pow(-1,lb-mb)*Angular::threej_2(2*lb,2*L,2*la,-twomb,twoM,twoma)*reduced_ME;
}

// Angular integral of the form 
// C(L,1,J;M-\lambda,\lambda,M) \Omega_b^T Y_{L,M-\lambda} sigma_i \Omega_a

// sigma_x
double ang_int_sig_x(int L, int twoJ, int twoM, int lambda, int lb, int twojb, int twomb, int la, int twoja, int twoma){

  double term_1 = cg(2*lb,1,twojb,twomb-1,1,twomb)*cg(2*la,1,twoja,twoma+1,-1,twoma)*ME_spherical_harm(L,twoM,lb,la,twomb-1,twoma+1);
  double term_2 = cg(2*lb,1,twojb,twomb+1,-1,twomb)*cg(2*la,1,twoja,twoma-1,1,twoma)*ME_spherical_harm(L,twoM,lb,la,twomb+1,twoma-1);

  return cg(2*L,2,twoJ,twoM-(2*lambda),2*lambda,twoM)*(term_1 + term_2);
}

// sigma_y
std::complex<double> ang_int_sig_y(int L, int twoJ, int twoM, int lambda, int lb, int twojb, int twomb, int la, int twoja, int twoma){

  double term_1 = -cg(2*lb,1,twojb,twomb-1,1,twomb)*cg(2*la,1,twoja,twoma+1,-1,twoma)*ME_spherical_harm(L,twoM,lb,la,twomb-1,twoma+1);
  double term_2 = cg(2*lb,1,twojb,twomb+1,-1,twomb)*cg(2*la,1,twoja,twoma-1,1,twoma)*ME_spherical_harm(L,twoM,lb,la,twomb+1,ma-1);
  
  return i*cg(2*L,2,twoJ,twoM-(2*lambda),2*lambda,twoM)*(term_1 + term_2);
}

// sigma_z
double ang_int_sig_z(int L, int twoJ, int twoM, int lambda, int lb, int twojb, int twomb, int la, int twoja, int twoma){

  double term_1 = cg(2*lb,1,twojb,twomb-1,1,twomb)*cg(2*la,1,twoja,twoma-1,1,twoma)*ME_spherical_harm(L,twoM,lb,la,twomb-1,twoma-1);
  double term_2 = -cg(2*lb,1,twojb,twomb+1,-1,twomb)*cg(2*la,1,twoja,twoma+1,-1,twoma)*ME_spherical_harm(L,twoM,lb,la,twomb+1,twoma+1);
  
  return cg(2*L,2,twoJ,twoM-(2*lambda),2*lambda,twoM)*(term_1 + term_2);
}

// Need to calculate the angular integral for three different vector spherical harmonic components
// Which will then be used to calculate the sum over sigma = -1, 0, 1 components

// \vec{Y}_{J,J,M}, L = J
// \vec{Y}_{J,J-1,M}, L = J - 1
// \vec{Y}_{J,J+1,M}, L = J + 1

// General function that calculates the entire angular and radial components for a vector spherical harmonic
// < b | alpha . jL \vec{Y}_{J,L,M} | a >
// for a given ma and mb (this will be used to calculate the reduced matrix element)

std::complex<double> ME_JLM(const DiracSpinor &Fa, const DiracSpinor &Fb, const auto jL, int L, int twoJ, int twoM, int twomb, int twoma){

  // Defining the quantum numbers to feed into the angular integrals
  auto la = Fa.l(); // int
  auto twoja = Fa.twoj(); // half-int
  auto l_til_a = Angular::l_tilde_k(Fa.kappa()); // int

  auto lb = Fb.l(); // int
  auto twojb = Fb.twoj(); // half-int
  auto l_til_b = Angular::l_tilde_k(Fb.kappa());

  const auto &grid = Fa.grid();
  using namespace qip::overloads;

  double radial_int_fbga = radial_int(Fb.f(), Fa.g(), jL, grid);
  double radial_int_gbfa = radial_int(Fb.g(), Fa.f(), jL, grid);

  std::complex<double> term1 = (1/sqrt(2))*i*radial_int_fbga*ang_int_sig_x(L,twoJ,twoM,-1,lb,twojb,twomb,l_til_a,twoja,twoma);
  std::complex<double> term2 = (1/sqrt(2))*radial_int_fbga*ang_int_sig_y(L,twoJ,twoM,-1,lb,twojb,twomb,l_til_a,twoja,twoma);
  std::complex<double> term3 = radial_int_fbga*ang_int_sig_z(L,twoJ,twoM,0,lb,twojb,twomb,l_til_a,twoja,twoma);
  std::complex<double> term4 = -(1/sqrt(2))*i*radial_int_fbga*ang_int_sig_z(L,twoJ,twoM,1,lb,twojb,twomb,l_til_a,twoja,twoma);
  std::complex<double> term5 = (1/sqrt(2))*radial_int_fbga*ang_int_sig_y(L,twoJ,twoM,1,lb,twojb,twomb,l_til_a,twoja,twoma);

  std::complex<double> term6 = -(1/sqrt(2))*i*radial_int_gbfa*ang_int_sig_x(L,twoJ,twoM,-1,l_til_b,twojb,twomb,la,twoja,twoma);
  std::complex<double> term7 = -(1/sqrt(2))*radial_int_gbfa*ang_int_sig_y(L,twoJ,twoM,-1,l_til_b,twojb,twomb,la,twoja,twoma);
  std::complex<double> term8 = -i*radial_int_gbfa*ang_int_sig_z(L,twoJ,twoM,0,l_til_b,twojb,twomb,la,twoja,twoma);
  std::complex<double> term9 = (1/sqrt(2))*radial_int_gbfa*ang_int_sig_x(L,twoJ,twoM,1,l_til_b,twojb,twomb,la,twoja,twoma);
  std::complex<double> term10 = -(1/sqrt(2))*radial_int_gbfa*ang_int_sig_y(L,twoJ,twoM,1,l_til_b,twojb,twomb,la,twoja,twoma);

  return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10;
}
