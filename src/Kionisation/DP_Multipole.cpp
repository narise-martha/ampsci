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
double cg(double la, double lb, int L, int ma, int mb, int M){
  return Angular::cg_2(int(2*la), 2*ma, int(2*lb), 2*mb, 2*L, 2*M);
}

// Matrix elements of a spherical harmonic, Y_{L,M}
// < lb, mb | Y_{L,M} | la, ma >
double ME_spherical_harm(int L, int M, int lb, int la, int mb, int ma) {
    
    // Reduced Matrix element
    // < lb || C^L_M || la >
    reduced_ME = pow(-1,lb)*sqrt((2*lb+1)*(2*la+1))*Angular::threej_2(2*lb,2*L,2*la,0,0,0);

    // Y_{L,M} = sqrt((2L+1)/4pi) C^L_M
    return sqrt((2*L+1)/4*M_PI)*pow(-1,lb-mb)*Angular::threej_2(2*lb,2*L,2*la,-2*mb,2*M,2*ma)*reduced_ME;
}

// Angular integral of the form 
// C(L,1,J;M-\lambda,\lambda,M) \Omega_b^T Y_{L,M-\lambda} sigma_i \Omega_a

// sigma_x
double ang_int_sig_x(int L, int J, int M, int lambda, int lb, int jb, int mb, int la, int ja, int ma){

  double term_1 = cg(lb,0.5,jb,mb-0.5,0.5,mb)*cg(la,0.5,ja,ma+0.5,-0.5,ma)*ME_spherical_harm(L,M,lb,la,mb-0.5,ma+0.5);
  double term_2 = cg(lb,0.5,jb,mb+0.5,-0.5,mb)*cg(la,0.5,ja,ma-0.5,0.5,ma)*ME_spherical_harm(L,M,lb,la,mb+0.5,ma-0.5);

  return cg(L,1,J,M-lambda,lambda,M)*(term_1 + term_2);
}

// sigma_y
std::complex<double> ang_int_sig_y(int L, int J, int M, int lambda, int lb, int jb, int mb, int la, int ja, int ma){

  double term_1 = -cg(lb,0.5,jb,mb-0.5,0.5,mb)*cg(la,0.5,ja,ma+0.5,-0.5,ma)*ME_spherical_harm(L,M,lb,la,mb-0.5,ma+0.5);
  double term_2 = cg(lb,0.5,jb,mb+0.5,-0.5,mb)*cg(la,0.5,ja,ma-0.5,0.5,ma)*ME_spherical_harm(L,M,lb,la,mb+0.5,ma-0.5);
  
  return i*cg(L,1,J,M-lambda,lambda,M)*(term_1 + term_2);
}

// sigma_z
double ang_int_sig_z(int L, int J, int M, int lambda, int lb, int jb, int mb, int la, int ja, int ma){

  double term_1 = cg(lb,0.5,jb,mb-0.5,0.5,mb)*cg(la,0.5,ja,ma-0.5,0.5,ma)*ME_spherical_harm(L,M,lb,la,mb-0.5,ma-0.5);
  double term_2 = -cg(lb,0.5,jb,mb+0.5,-0.5,mb)*cg(la,0.5,ja,ma+0.5,-0.5,ma)*ME_spherical_harm(L,M,lb,la,mb+0.5,ma+0.5);
  
  return cg(L,1,J,M-lambda,lambda,M)*(term_1 + term_2);
}

// Need to calculate the angular integral for three different vector spherical harmonic components
// Which will then be used to calculate the sum over sigma = -1, 0, 1 components

// \vec{Y}_{J,J,M}, L = J
// \vec{Y}_{J,J-1,M}, L = J - 1
// \vec{Y}_{J,J+1,M}, L = J + 1

// General function that calculates the entire angular and radial components and then sums them up
// < b || alpha . jL \vec{Y}_{J,L,M} || a >
// for a given ma and mb (this will be used to calculate the reduced matrix element)
std::complex<double> vec_sph_LM(const DiracSpinor &Fa, const DiracSpinor &Fb, const auto jL, int L, int J, int M, int mb, int la, int ja, int ma){

    // Defining the quantum numbers to feed into the angular integrals
  auto la = Fa.l();
  auto twoja = Fa.twoj();
  auto l_til_a = Angular::l_tilde_k(Fa.kappa());

  auto lb = Fb.l();
  auto twojb = Fb.twoj();
  auto l_til_b = Angular::l_tilde_k(Fb.kappa());

  double term1 = (1/sqrt(2))*i*radial_int()*ang_int_sig_x(L,J,M,-1,lb,jb,mb,la,ja,ma)

  return ;
}
