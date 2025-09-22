#include "Angular/Wigner369j.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/DiagramRPA0_jL.hpp"
#include "HF/HartreeFock.hpp"
#include "Kion_functions.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/UnitConv_conversions.hpp"
// #include "StandardHaloModel.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/ostream.hpp"
#include "qip/Methods.hpp"
#include "qip/omp.hpp"
#include <iostream>
#include <memory>
#include <complex>
#include <optional>

//==============================================================================

// RADIAL INTEGRAL FUNCTION
// int func1 * jL * func2 dr
double radial_int(const std::vector<double> &func1,
                  const std::vector<double> &func2,
                  const std::vector<double> &jL, Grid grid) {
  return NumCalc::integrate(grid.du(), 0, 0, grid.drdu(), func1, func2, jL);
}

//==============================================================================

// ANGULAR INTAGRAL FUNCTIONS

// Clecbsch-Gordon coefficient but in a different order
double cg(const int twola, const int twolb, const int twoL, const int twoma, const int twomb, const int twoM){
  return Angular::cg_2(twola, twoma, twolb, twomb, twoL, twoM);
}

// Matrix elements of a spherical harmonic, Y_{L,M}
// < lb, mb | Y_{L,M} | la, ma >
double ME_spherical_harm(const int twoL, const int twoM, const int lb, const int la, const int twomb, const int twoma) {
    
    // int J = int(0.5*( twoJ + 0.001));

    // Reduced Matrix element
    // < lb || C^L_M || la >
    double reduced_ME = pow(-1,lb)*sqrt((2*lb+1)*(2*la+1))*Angular::threej_2(2*lb,twoL,2*la,0,0,0);
    
    // (sqrt((twoJ+1)/(4*M_PI)))*Angular::Ck_kk(J,-kappa_b,kappa_a); 
    
    //pow(-1,lb)*sqrt((2*lb+1)*(2*la+1))*Angular::threej_2(2*lb,twoL,2*la,0,0,0);

    // Y_{L,M} = sqrt((2L+1)/4pi) C^L_M
    // pow(-1,0.5*(2*lb-(twomb+0.001)))*
    // auto sign = Angular::neg1pow_2(2*lb - twomb)
    return sqrt((twoL+1)/4*M_PI)*Angular::neg1pow_2(2*lb - twomb)*Angular::threej_2(2*lb,twoL,2*la,-twomb,twoM,twoma)*reduced_ME;
}

// Angular integral of the form 
// C(L,1,J;M-\lambda,\lambda,M) \Omega_b^T Y_{L,M-\lambda} sigma_i \Omega_a

// sigma_x
double ang_int_sig_x(const int twoL, const int twoJ, const int twoM, const int lambda, const int lb, const int twojb, const int twomb, const int la, const int twoja, const int twoma){

  double term_1 = cg(2*lb,1,twojb,twomb-1,1,twomb)*cg(2*la,1,twoja,twoma+1,-1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb-1,twoma+1);
  double term_2 = cg(2*lb,1,twojb,twomb+1,-1,twomb)*cg(2*la,1,twoja,twoma-1,1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb+1,twoma-1);

  return cg(twoL,2,twoJ,twoM-(2*lambda),2*lambda,twoM)*(term_1 + term_2);
}

// sigma_y
std::complex<double> ang_int_sig_y(const int twoL, const int twoJ, const int twoM, const int lambda, const int lb, const int twojb, const int twomb, const int la, const int twoja, const int twoma){

  std::complex<double> i(0.0,1.0); 

  double term_1 = -cg(2*lb,1,twojb,twomb-1,1,twomb)*cg(2*la,1,twoja,twoma+1,-1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb-1,twoma+1);
  double term_2 = cg(2*lb,1,twojb,twomb+1,-1,twomb)*cg(2*la,1,twoja,twoma-1,1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb+1,twoma-1);
  
  return i*cg(twoL,2,twoJ,twoM-(2*lambda),2*lambda,twoM)*(term_1 + term_2);
}

// sigma_z
double ang_int_sig_z(const int twoL, const int twoJ, const int twoM, const int lambda, const int lb, const int twojb, const int twomb, const int la, const int twoja, const int twoma){

  double term_1 = cg(2*lb,1,twojb,twomb-1,1,twomb)*cg(2*la,1,twoja,twoma-1,1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb-1,twoma-1);
  double term_2 = -cg(2*lb,1,twojb,twomb+1,-1,twomb)*cg(2*la,1,twoja,twoma+1,-1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb+1,twoma+1);
  
  return cg(twoL,2,twoJ,twoM-(2*lambda),2*lambda,twoM)*(term_1 + term_2);
}

// Need to calculate the angular integral for three different vector spherical harmonic components
// Which will then be used to calculate the sum over sigma = -1, 0, 1 components

// \vec{Y}_{J,J,M}, L = J
// \vec{Y}_{J,J-1,M}, L = J - 1
// \vec{Y}_{J,J+1,M}, L = J + 1

// General function that calculates the entire angular and radial components for a vector spherical harmonic
// < b | alpha . jL \vec{Y}_{J,L,M} | a >
// for a given ma and mb (this will be used to calculate the reduced matrix element)

std::complex<double> ME_jL_Y_JLM(const DiracSpinor &Fa, const DiracSpinor &Fb, const std::vector<double> jL, const int twoL, const int twoJ, const int twoM, const int twomb, const int twoma){

  // Defining the quantum numbers to feed into the angular integrals
  auto la = Fa.l(); // int
  auto twoja = Fa.twoj(); // half-int
  auto l_til_a = Angular::l_tilde_k(Fa.kappa()); // int

  auto lb = Fb.l(); // int
  auto twojb = Fb.twoj(); // half-int
  auto l_til_b = Angular::l_tilde_k(Fb.kappa()); // int

  const auto &grid = Fa.grid();
  using namespace qip::overloads;

  double radial_int_fbga = radial_int(Fb.f(), Fa.g(), jL, grid);
  double radial_int_gbfa = radial_int(Fb.g(), Fa.f(), jL, grid);

  std::complex<double> i(0.0,1.0); 

  std::complex<double> term1 = (1/sqrt(2))*i* radial_int_fbga*ang_int_sig_x(twoL,twoJ,twoM,-1,lb,twojb,twomb,l_til_a,twoja,twoma);
  std::complex<double> term2 = (1/sqrt(2))*   radial_int_fbga*ang_int_sig_y(twoL,twoJ,twoM,-1,lb,twojb,twomb,l_til_a,twoja,twoma);
  std::complex<double> term3 =                radial_int_fbga*ang_int_sig_z(twoL,twoJ,twoM, 0,lb,twojb,twomb,l_til_a,twoja,twoma);
  std::complex<double> term4 = -(1/sqrt(2))*i*radial_int_fbga*ang_int_sig_z(twoL,twoJ,twoM, 1,lb,twojb,twomb,l_til_a,twoja,twoma);
  std::complex<double> term5 = (1/sqrt(2))*   radial_int_fbga*ang_int_sig_y(twoL,twoJ,twoM, 1,lb,twojb,twomb,l_til_a,twoja,twoma);

  std::complex<double> term6 = -(1/sqrt(2))*i*radial_int_gbfa*ang_int_sig_x(twoL,twoJ,twoM,-1,l_til_b,twojb,twomb,la,twoja,twoma);
  std::complex<double> term7 = -(1/sqrt(2))*  radial_int_gbfa*ang_int_sig_y(twoL,twoJ,twoM,-1,l_til_b,twojb,twomb,la,twoja,twoma);
  std::complex<double> term8 = -i*            radial_int_gbfa*ang_int_sig_z(twoL,twoJ,twoM, 0,l_til_b,twojb,twomb,la,twoja,twoma);
  std::complex<double> term9 = (1/sqrt(2))*   radial_int_gbfa*ang_int_sig_x(twoL,twoJ,twoM, 1,l_til_b,twojb,twomb,la,twoja,twoma);
  std::complex<double> term10 = -(1/sqrt(2))* radial_int_gbfa*ang_int_sig_y(twoL,twoJ,twoM, 1,l_til_b,twojb,twomb,la,twoja,twoma);

  return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10;
}

//==============================================================================

// Calculating |< jb || alpha . a^(sigma)_J || ja >|^2
// but using ME calculated by Johnson. 

// Calculating <ka||C^k||kb> using Angular::Ck_kk(int k, int ka, int kb)
// Rescaling from C^sigma -> Y^sigma
double RME_sigma_Y(int kappa_b, int kappa_a, const int sigma, const int twoJ){
  
  int J = int(0.5*( twoJ + 0.001)); // Convert int twoJ -> int J

  // Matrix elements calculated by Johnson (p. 143)
  if (sigma == -1){
    return -(sqrt((twoJ+1)/(4*M_PI)))*Angular::Ck_kk(J,-kappa_b,kappa_a);
  }

  else if (sigma == 0){
    // Ensuring no divide by zero errors occur
    if (twoJ==-1 || twoJ == 0){
      return 0;
    }

    else{
    return (sqrt((twoJ+1)/(4*M_PI)))
          *((kappa_a - kappa_b)/(sqrt(J*(J+1))))*Angular::Ck_kk(J,kappa_b,kappa_a);
    }
  }

  else if (sigma == 1){
    // Ensuring no divide by zero errors occur
    if (twoJ==-1 || twoJ == 0){
      return 0;
    }
    else{
      return (sqrt((twoJ+1)/(4*M_PI)))*((kappa_a + kappa_b)/(sqrt(J*(J+1))))*Angular::Ck_kk(J,-kappa_b,kappa_a);
    }
  }
}

std::complex<double> RME_alpha_jL_Y(const DiracSpinor &Fa, const DiracSpinor &Fb, const int sigma, const int twoJ, const std::vector<double> jL){

  const auto &grid = Fa.grid();
  using namespace qip::overloads;

  std::complex<double> i(0.0,1.0); 

  return i*radial_int(Fb.f(), Fa.g(), jL, grid)*RME_sigma_Y(Fb.kappa(), -Fa.kappa(),sigma,twoJ) 
        - i*radial_int(Fb.g(), Fa.f(), jL, grid)*RME_sigma_Y(-Fb.kappa(), Fa.kappa(),sigma,twoJ);
}






// < b | alpha . a^(sigma)_{J,M} | a >
// Where sigma = -1,0,1
// a^(-1)_{J,M} = sqrt(J/(2J+1)) a_{J,J-1,M} - sqrt((J+1)/(2J+1)) a_{J,J+1,M}
// a^(0)_{J,M} = a_{J,J,M}
// a^(1)_{J,M} = sqrt((J+1)/(2J+1)) a_{J,J-1,M} + sqrt(J/(2J+1)) a_{J,J+1,M}

// Using the Wigner-Eckhart theorem, 
// < b | alpha . a^(sigma)_{J,M} | a > = (-1)^(jb - ja) Wig3j(jb, M, ja; -mb, J, ma) <jb || alpha . a^(sigma)_J || ja >

//==============================================================================

// Reduced matrix element
// < jb || alpha . a^(sigma)_J || ja > = (-1)^(ja-jb) [Wig3j(jb, J, ja; -1/2, 0, 1/2)]^{-1} < jb, 1/2 | alpha . a^(sigma)_{J,0} | ja, 1/2 >
// absolute value squared => |< jb || alpha . a^(sigma)_J || ja >|^2
// For a single state (need to sum over initial and final states)
double abs_RME_alpha_a(bool Johnson, const DiracSpinor &Fa, const DiracSpinor &Fb, const std::vector<double> jJ, const std::vector<double> jJplus1, const std::vector<double> jJminus1, const int sigma, const int twoJ){

  double eta = double(sqrt(twoJ/(2*twoJ+2)));
  double zeta = double(sqrt((twoJ+2)/(2*twoJ+2)));

  // Using my matrix elements
  if (Johnson == false){

    std::complex<double> full_ME;

    if (sigma == -1){
      // a^(-1)_{J,M} = - eta a_{J,J-1,M} - zeta a_{J,J+1,M}
      full_ME = -eta*ME_jL_Y_JLM(Fa,Fb,jJminus1,twoJ-2,twoJ,0,1,1) // L = J - 1
              - zeta*ME_jL_Y_JLM(Fa,Fb,jJplus1,twoJ+2,twoJ,0,1,1); // L = J + 1
    }

    else if (sigma == 0){
      // a^(0)_{J,M} = a_{J,J,M}
      full_ME = ME_jL_Y_JLM(Fa,Fb,jJ,twoJ,twoJ,0,1,1);  // L = J
    }

    else if (sigma == 1){
      // a^(1)_{J,M} = zeta a_{J,J-1,M} - eta a_{J,J+1,M}

      full_ME = zeta*ME_jL_Y_JLM(Fa,Fb,jJminus1,twoJ-2,twoJ,0,1,1) // L = J - 1
              - eta*ME_jL_Y_JLM(Fa,Fb,jJplus1,twoJ+2,twoJ,0,1,1); // L = J + 1
    }

    if (Angular::threej_2(Fb.twoj(),twoJ,Fa.twoj(),-1,0,1) == 0){
    return 0; 
  }

    // Using Wigner-Ekhart to convert from full matrix element to reduced matrix element
    return pow(abs(Angular::neg1pow_2(1 - Fb.twoj())
            *(1/Angular::threej_2(Fb.twoj(),twoJ,Fa.twoj(),-1,0,1)) * full_ME),2);
  }

  // USING JOHNSON matrix elements
  // Instead of calculating the full matrix element, we already
  // know the reduced matrix elements (since they can be composed of
  //  < b || C_J || a >, which are just rescaled spherical harmonics)

  if (Johnson == true){

    std::complex<double> full_RME;

    if (sigma == -1){

      // a^(-1)_{J,M} = sqrt(J/(2J+1)) j_(J-1) Y_{J,J-1,M} - sqrt((J+1)/(2J+1)) j_(J+1) Y_{J,J+1,M}

      full_RME = -eta*eta*RME_alpha_jL_Y(Fa,Fb,-1,twoJ,jJminus1) 
                - eta*zeta*RME_alpha_jL_Y(Fa,Fb,1,twoJ,jJminus1) 
                + zeta*zeta*RME_alpha_jL_Y(Fa,Fb,-1,twoJ,jJplus1) 
                - eta*zeta*RME_alpha_jL_Y(Fa,Fb,1,twoJ,jJplus1);
    }

    else if (sigma == 0){

      // a^(0)_{J,M} = jJ Y_{J,J,M}

      full_RME = RME_alpha_jL_Y(Fa,Fb,0,twoJ,jJ);
    }

    else if (sigma == -1){

      // a^(1)_{J,M} = sqrt((J+1)/(2J+1)) j_(J-1) Y_{J,J-1,M} + sqrt(J/(2J+1)) j_(J+1) Y_{J,J+1,M}

      full_RME = -eta*zeta*RME_alpha_jL_Y(Fa,Fb,-1,twoJ,jJminus1) 
                - zeta*zeta*RME_alpha_jL_Y(Fa,Fb,1,twoJ,jJminus1) 
                - eta*zeta*RME_alpha_jL_Y(Fa,Fb,-1,twoJ,jJplus1) 
                + eta*eta*RME_alpha_jL_Y(Fa,Fb,1,twoJ,jJplus1);
    }

    return pow(abs(full_RME),2);
  } 

  }




//////////////////////////////////////////////////////////////////////////////////////////////
// This function calculates the reduced matrix element (RME) of
// < jb || alpha . a^(sigma)_J || ja >
// for a given initial state, and sums over all the final states that result in an ionised state.
// E = mc^2, energy of the dark photon, m, mass of dark photon 
double abs_RME_total(bool Johnson, const HF::HartreeFock *vHF, const DiracSpinor &Fa, const double E, const int twoJ, const std::vector<double> jJ, const std::vector<double> jJplus1, const std::vector<double> jJminus1) {

  ContinuumOrbitals cntm(vHF);

  // Calculating final energy state.
  // E_f = mc^2 - |E_i| = mc^2 + E_i
  // since E_i < 0
  const auto ec = E + Fa.en();

  // If final energy is negative, ignore this state
  // Only want ionised final states.
  if (ec < 0.0) {
    return 0.0;
  }

  // Solving for continuum states
  cntm.solveContinuumHF(ec, 0, 4, &Fa);

  // Calculating the RME of a single state, and then
  // summing over the final states states, and sigma
  double RME_tot = 0.0;
  for (const auto &Fb : cntm.orbitals) {
    for (int sigma = -1; sigma < 2; sigma++){
    RME_tot += abs_RME_alpha_a(Johnson,Fa, Fb, jJ, jJplus1, jJminus1, sigma, twoJ); 
    }
  }
  return RME_tot;
}

// Same as the function above, but we only want to sum over sigma = 0,1
// Could make this more compact by adding another Boolian in the above function
// For EDA true/false to change the sum over sigma
double abs_RME_total_EDA(bool Johnson, const HF::HartreeFock *vHF, const DiracSpinor &Fa, const double E, const int twoJ, const std::vector<double> jJ, const std::vector<double> jJplus1, const std::vector<double> jJminus1) {

  ContinuumOrbitals cntm(vHF);

  // Calculating final energy state.
  // E_f = mc^2 - |E_i| = mc^2 + E_i
  // since E_i < 0
  const auto ec = E + Fa.en();

  // If final energy is negative, ignore this state
  // Only want ionised final states.
  if (ec < 0.0) {
    return 0.0;
  }

  // Solving for continuum states
  cntm.solveContinuumHF(ec, 0, 4, &Fa);

  // Calculating the RME of a single state, and then
  // summing over the final states states, and sigma
  double RME_tot = 0.0;
  for (const auto &Fb : cntm.orbitals) {
    for (int sigma = 0; sigma < 2; sigma++){
    RME_tot += abs_RME_alpha_a(Johnson,Fa, Fb, jJ, jJplus1, jJminus1, sigma, twoJ); 
    }
  }
  return RME_tot;
}

/*
// Taking the sum over sigma for the RME, and multiplying by 4pi
double abs_Tab_sq(const DiracSpinor &Fa, const DiracSpinor &Fb){

  double summed_reduced_ME = 0;
  
  for (int sigma = -1; sigma < 2; sigma++){
    for (int twoJ = ; twoJ < ; twoJ++){
      summed_reduced_ME += abs_RME_total(vHF, Fa, E, L, jL)
    }
  }

  return 4*M_PI*summed_reduced_ME
}

 
std::complex<double> ME_WignerEkhart(){

  std::complex<double> reduced_ME = (1/Angular::threej_2(twojb,twoM,twoja,-1,twoJ,1))*ME_Y_JLM(Fa,Fb,jL,twoL,twoJ,twoM,1,1);

  return Angular::threej_2(twojb,twoM,twoja,-twoma,twoJ,)
}

// I think all of this is kind of wrong?
std::complex<double> ME_a_JM(const DiracSpinor &Fa, const DiracSpinor &Fb, int sigma, int twoJ, int twoM){

  if (sigma==-1){
    return sqrt(twoJ/(2*twoJ+2))*ME_Y_JLM(Fa,Fb,jL,twoJ-2,twoJ,twoM,1,1) - sqrt((twoJ+1)/(2*twoJ+2))*ME_Y_JLM(Fa,Fb,jL,twoJ+2,twoJ,twoM,1,1);
  }

  elseif (sigma == 0){
    return ME_Y_JLM(Fa,Fb,jL,twoJ,twoJ,twoM,1,1);
  }

  elseif (sigma == 1){
    return sqrt((twoJ+1)/(2*twoJ+2))*ME_Y_JLM(Fa,Fb,jL,twoJ-2,twoJ,twoM,1,1) + sqrt(twoJ/(2*twoJ+2))*ME_Y_JLM(Fa,Fb,jL,twoJ+2,twoJ,twoM,1,1);
  }

  else{
    return std::nullopt;
  }
}

*/
