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

double R_plus(const DiracSpinor &Fa, const DiracSpinor &Fb, const std::vector<double> &tK, Grid grid){
  return radial_int(Fa.f(), Fb.f(), tK, grid) + radial_int(Fa.g(), Fb.g(), tK, grid);
}

double R_minus(const DiracSpinor &Fa, const DiracSpinor &Fb, const std::vector<double> &tK, Grid grid){
  return radial_int(Fa.f(), Fb.f(), tK, grid) - radial_int(Fa.g(), Fb.g(), tK, grid);
}

double P_plus(const DiracSpinor &Fa, const DiracSpinor &Fb, const std::vector<double> &tK, Grid grid){
  return radial_int(Fb.f(), Fa.g(), tK, grid) + radial_int(Fb.g(), Fa.f(), tK, grid);
}

double P_minus(const DiracSpinor &Fa, const DiracSpinor &Fb, const std::vector<double> &tK, Grid grid){
  return radial_int(Fa.f(), Fb.g(), tK, grid) - radial_int(Fa.g(), Fb.f(), tK, grid);
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
    return sqrt((twoL+1)/4*M_PI)*Angular::neg1pow_2(2*lb - twomb)*Angular::threej_2(2*lb,twoL,2*la,-twomb,twoM,twoma)*reduced_ME;
}

//==============================================================================

// My version of calculations

// Angular integral of the form 
// C(L,1,J;M-\lambda,\lambda,M) \Omega_b^T Y_{L,M-\lambda} sigma_i \Omega_a
double ang_int(const char i, const int twoL, const int twoJ, const int lambda, const int lb, const int twojb, const int la, const int twoja){

  int twoma = 1;
  int twomb = 1;
  int twoM = 0;
  double term_1;
  double term_2;

  if (i == 'x'){
    term_1 = cg(2*lb,1,twojb,twomb-1,1,twomb)*cg(2*la,1,twoja,twoma+1,-1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb-1,twoma+1);
    term_2 = cg(2*lb,1,twojb,twomb+1,-1,twomb)*cg(2*la,1,twoja,twoma-1,1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb+1,twoma-1);
  }

  else if (i == 'y'){
    // This technically is imaginary, but i is taken out and put in the sigma_y terms ME_jL_Y_JLM
    term_1 = -cg(2*lb,1,twojb,twomb-1,1,twomb)*cg(2*la,1,twoja,twoma+1,-1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb-1,twoma+1);
    term_2 = cg(2*lb,1,twojb,twomb+1,-1,twomb)*cg(2*la,1,twoja,twoma-1,1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb+1,twoma-1);
  }

  else if (i == 'z'){
    term_1 = cg(2*lb,1,twojb,twomb-1,1,twomb)*cg(2*la,1,twoja,twoma-1,1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb-1,twoma-1);
    term_2 = -cg(2*lb,1,twojb,twomb+1,-1,twomb)*cg(2*la,1,twoja,twoma+1,-1,twoma)*ME_spherical_harm(twoL,twoM-2*lambda,lb,la,twomb+1,twoma+1);
  }

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

std::complex<double> ME_jL_Y_JLM(const DiracSpinor &Fa, const DiracSpinor &Fb, const std::vector<double> jL, const int twoL, const int twoJ){

  // Defining the quantum numbers to feed into the angular integrals
  auto la = Fa.l(); // int
  auto twoja = Fa.twoj(); // j is half-int
  auto l_til_a = Angular::l_tilde_k(Fa.kappa()); // int

  auto lb = Fb.l(); // int
  auto twojb = Fb.twoj(); // j is half-int
  auto l_til_b = Angular::l_tilde_k(Fb.kappa()); // int

  const auto &grid = Fa.grid();
  using namespace qip::overloads;

  double radial_int_fbga = radial_int(Fb.f(), Fa.g(), jL, grid);  // V_{ij} [jL]
  double radial_int_gbfa = radial_int(Fb.g(), Fa.f(), jL, grid);  // W_{ij} [jL]

  std::complex<double> i(0.0,1.0); 

  // There are so many terms because of the sum over the polarisation vectors, dot producted with the
  // spin-1 basis vactors. 
  std::complex<double> term1 = (1/sqrt(2))*i* radial_int_fbga*ang_int('x',twoL,twoJ,-1,lb,twojb,l_til_a,twoja);
  std::complex<double> term2 = (1/sqrt(2))*i* radial_int_fbga*ang_int('y',twoL,twoJ,-1,lb,twojb,l_til_a,twoja);
  std::complex<double> term3 =                radial_int_fbga*ang_int('z',twoL,twoJ, 0,lb,twojb,l_til_a,twoja);
  std::complex<double> term4 = -(1/sqrt(2))*i*radial_int_fbga*ang_int('x',twoL,twoJ, 1,lb,twojb,l_til_a,twoja);
  std::complex<double> term5 = (1/sqrt(2))*i*  radial_int_fbga*ang_int('y',twoL,twoJ,1,lb,twojb,l_til_a,twoja);

  std::complex<double> term6 = -(1/sqrt(2))*i*radial_int_gbfa*ang_int('x',twoL,twoJ,-1,l_til_b,twojb,la,twoja);
  std::complex<double> term7 = -(1/sqrt(2))*i*radial_int_gbfa*ang_int('y',twoL,twoJ,-1,l_til_b,twojb,la,twoja);
  std::complex<double> term8 = -i*            radial_int_gbfa*ang_int('z',twoL,twoJ, 0,l_til_b,twojb,la,twoja);
  std::complex<double> term9 = (1/sqrt(2))*   radial_int_gbfa*ang_int('x',twoL,twoJ, 1,l_til_b,twojb,la,twoja);
  std::complex<double> term10 = -(1/sqrt(2))*i*radial_int_gbfa*ang_int('y',twoL,twoJ,1,l_til_b,twojb,la,twoja);

  return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10;
}

//==============================================================================

// JOHNSON METHOD

// Calculating |< jb || alpha . a^(sigma)_J || ja >|^2

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
    if (J==-1 || J == 0){
      return 0;
    }

    else{
    return (sqrt((twoJ+1)/(4*M_PI)))
          *((kappa_a - kappa_b)/(sqrt(J*(J+1))))*Angular::Ck_kk(J,kappa_b,kappa_a);
    }
  }

  else if (sigma == 1){
    // Ensuring no divide by zero errors occur
    if (J==-1 || J == 0){
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

//==============================================================================



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
double abs_RME_alpha_a(int method, const DiracSpinor &Fa, const DiracSpinor &Fb, 
                        const std::vector<double> jJ, const std::vector<double> jJplus1, 
                        const std::vector<double> jJminus1, const int sigma, const int twoJ, 
                        const std::vector<double> qr){

  double eta = double(sqrt(twoJ/(2*twoJ+2)));
  double zeta = double(sqrt((twoJ+2)/(2*twoJ+2)));

  int J = int(0.5*( twoJ + 0.001)); // Convert int twoJ -> int J
  const auto &grid = Fa.grid();
  using namespace qip::overloads;
  const auto &gr = Fb.grid();

  std::complex<double> i(0.0,1.0);
  std::vector<double> jJ_dash = J*(jJ/qr) - jJplus1; //NumCalc::derivative(jJ,gr.drdu(), gr.du(), 1);

  double temporal_RME = R_plus(Fa,Fb,jJ,grid)*Angular::Ck_kk(J,Fb.kappa(),Fa.kappa());

  // Using my matrix elements
  if (method == 0){

    std::complex<double> spatial_ME;

    if (sigma == -1){
      // a^(-1)_{J,M} = - eta a_{J,J-1,M} - zeta a_{J,J+1,M}
      spatial_ME = eta*ME_jL_Y_JLM(Fa,Fb,jJminus1,twoJ-2,twoJ) // L = J - 1
              - zeta*ME_jL_Y_JLM(Fa,Fb,jJplus1,twoJ+2,twoJ); // L = J + 1
    }

    else if (sigma == 0){
      // a^(0)_{J,M} = a_{J,J,M}
      spatial_ME = ME_jL_Y_JLM(Fa,Fb,jJ,twoJ,twoJ);  // L = J
    }

    else if (sigma == 1){
      // a^(1)_{J,M} = zeta a_{J,J-1,M} - eta a_{J,J+1,M}

      spatial_ME = zeta*ME_jL_Y_JLM(Fa,Fb,jJminus1,twoJ-2,twoJ) // L = J - 1
              + eta*ME_jL_Y_JLM(Fa,Fb,jJplus1,twoJ+2,twoJ); // L = J + 1
    }

    if (Angular::threej_2(Fb.twoj(),twoJ,Fa.twoj(),-1,0,1) == 0){
    return 0; 
  }

    // Using Wigner-Ekhart to convert from full matrix element to reduced matrix element
    return pow(-abs(Angular::neg1pow_2(1 - Fb.twoj())
            *(1/Angular::threej_2(Fb.twoj(),twoJ,Fa.twoj(),-1,0,1)) * spatial_ME + temporal_RME),2);
  }

  // USING JOHNSON matrix elements
  // Instead of calculating the full matrix element, we already
  // know the reduced matrix elements (since they can be composed of
  //  < b || C_J || a >, which are just rescaled spherical harmonics)

  if (method == 1){

    std::complex<double> full_RME;

    if (sigma == -1){

      // a^(-1)_{J,M} = eta j_(J-1) Y_{J,J-1,M} - zeta j_(J+1) Y_{J,J+1,M}

      full_RME = -eta*eta*RME_alpha_jL_Y(Fa,Fb,-1,twoJ,jJminus1) 
                - eta*zeta*RME_alpha_jL_Y(Fa,Fb,1,twoJ,jJminus1) 
                + zeta*zeta*RME_alpha_jL_Y(Fa,Fb,-1,twoJ,jJplus1) 
                - eta*zeta*RME_alpha_jL_Y(Fa,Fb,1,twoJ,jJplus1);
    }

    else if (sigma == 0){

      // a^(0)_{J,M} = jJ Y_{J,J,M}

      full_RME = RME_alpha_jL_Y(Fa,Fb,0,twoJ,jJ);
    }

    else if (sigma == 1){

      // a^(1)_{J,M} = zeta j_(J-1) Y_{J,J-1,M} + eta j_(J+1) Y_{J,J+1,M}

      full_RME = -eta*zeta*RME_alpha_jL_Y(Fa,Fb,-1,twoJ,jJminus1) 
                - zeta*zeta*RME_alpha_jL_Y(Fa,Fb,1,twoJ,jJminus1) 
                - eta*zeta*RME_alpha_jL_Y(Fa,Fb,-1,twoJ,jJplus1) 
                + eta*eta*RME_alpha_jL_Y(Fa,Fb,1,twoJ,jJplus1);
    }

    return pow(abs(-full_RME + temporal_RME),2);
  } 

  // Using Ben's calculations from the overleaf
  // For vector case, gamma^\tila = 1

  if (method == 2){

    std::complex<double> spatial_RME;

    if (sigma == -1){
      spatial_RME = i*(P_minus(Fa,Fb,jJ_dash,grid)-(Fb.kappa() - Fa.kappa())
                  *P_plus(Fa,Fb,jJ/qr,grid))*Angular::Ck_kk(J,Fb.kappa(),Fa.kappa());
    }

    if (sigma == 0){

      if (J==-1 || J == 0){
      return 0;
      }

      spatial_RME = -i*P_plus(Fa,Fb,jJ,grid)*((Fb.kappa() + Fa.kappa())/sqrt(J*(J+1)))
                  *Angular::Ck_kk(J,Fb.kappa(),-Fa.kappa());
    }

    if (sigma == 1){
      
      if (J==-1 || J == 0){
      return 0;
      }

      spatial_RME = i*(sqrt(J*(J+1))*P_minus(Fa,Fb,jJ/qr,grid) - ((Fb.kappa() - Fa.kappa())/sqrt(J*(J+1)))
                  *P_plus(Fa,Fb,(jJ/qr)+ jJ_dash,grid))*Angular::Ck_kk(J,Fb.kappa(),Fa.kappa());
    }

    return (2*J + 1)*pow(abs(-spatial_RME + 0*temporal_RME),2);
  }

  // Using the equations directly from Johnson (p. 190) for the transverse (velocity) gauge

  if (method == 3){
    std::complex<double> total_RME;

    if (sigma == -1){
      total_RME = 0;
    }

    if (sigma == 0){
      
      if (J==-1){
      return 0;
      }

      total_RME = ((Fb.kappa() + Fa.kappa())/(J+1))*Angular::Ck_kk(J,-Fb.kappa(),Fa.kappa())
                  *P_plus(Fa,Fb,jJ,grid);
    }

    if (sigma == 1){
      
      if (J==-1){
      return 0;
      }

      total_RME = Angular::Ck_kk(J,Fb.kappa(),Fa.kappa())
                *(((Fa.kappa()-Fb.kappa())/(J+1))*P_plus(Fa,Fb,jJ_dash+(jJ/qr),grid) 
                - J*P_minus(Fa,Fb,jJ/qr,grid)); // I made this a minus, because I think Johnson's integrals are the other way around
    }

    if (J == 0){
      return 0;
      }

    return pow(abs(i*sqrt((2*J + 1)*(J+1)/(4*M_PI*J))*total_RME),2);

  }

  }

//==============================================================================

// This function calculates the reduced matrix element (RME) of
// < jb || alpha . a^(sigma)_J || ja >
// for a given initial state, and sums over all the final states that result in an ionised state.
// E = mc^2, energy of the dark photon, m, mass of dark photon 
double abs_RME_total(int method, const HF::HartreeFock *vHF, const DiracSpinor &Fa, 
                    const double E, const int twoJ, const std::vector<double> jJ, 
                    const std::vector<double> jJplus1, const std::vector<double> jJminus1, 
                    const std::vector<double> qr) {

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
      double omega = (Fb.en() - Fa.en());

      RME_tot += (1/pow(omega,1))*abs_RME_alpha_a(method,Fa, Fb, jJ, jJplus1, jJminus1, sigma, twoJ,qr); 
    }
  }
  return RME_tot;
}

// Same as the function above, but we only want to sum over sigma = 0,1
// Could make this more compact by adding another Boolian in the above function
// For EDA true/false to change the sum over sigma
double abs_RME_total_EDA(int method, const HF::HartreeFock *vHF, const DiracSpinor &Fa, 
                          const double E, const int twoJ, const std::vector<double> jJ, 
                          const std::vector<double> jJplus1, const std::vector<double> jJminus1, 
                          const std::vector<double> qr) {

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
      double omega = (Fb.en() - Fa.en()); 

      RME_tot += (1/pow(omega,1))*abs_RME_alpha_a(method,Fa, Fb, jJ, jJplus1, jJminus1, sigma, twoJ,qr); 
    }
  }
  return RME_tot;
}

double abs_RME_total_EDA_Ben(const HF::HartreeFock *vHF, const DiracSpinor &Fa, 
                          const double E, const std::vector<double> jJ) {

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
  const auto &grid = Fa.grid();
  using namespace qip::overloads;

  double RME;
  double temporal_RME;
  double omega;
  
  std::vector<double> ones(Fa.grid().size(),1);

  // Solving for continuum states
  cntm.solveContinuumHF(ec, 0, 4, &Fa);

  // Calculating the RME of a single state, and then
  // summing over the final states states, and sigma
  double RME_tot = 0.0;
  for (const auto &Fb : cntm.orbitals) {
      omega = (Fb.en() - Fa.en()); 
      RME = (sqrt(2)/3)*(P_minus(Fa,Fb,ones,grid) + (Fb.kappa() - Fa.kappa())*P_plus(Fa,Fb,ones,grid))*Angular::Ck_kk(1,Fb.kappa(),Fa.kappa());
      temporal_RME = R_plus(Fa,Fb,jJ,grid)*Angular::Ck_kk(1,Fb.kappa(),Fa.kappa());
      RME_tot += (1/omega)*pow(temporal_RME - RME,2);
  }
  return RME_tot;
}

//==============================================================================

// Electric dipole approximation
// LENGTH FORM
// K = |<b|z|a>|^2 = (1/3) omega^2 alpha^2 |<b||r||a>|^2
double RME_single_e(const DiracSpinor &Fa, const DiracSpinor &Fb, const auto r) {

  const auto &grid = Fa.grid();
  using namespace qip::overloads;

  // Radial component int (fa*fb*r + ga*gb*r) dr
  double R = radial_int(Fa.f(), Fb.f(), r, grid) + radial_int(Fa.g(), Fb.g(), r, grid);

  // Angular component <b|| C^L_q || a>
  // Electric dipole approximation, L=1
  double C = Angular::Ck_kk(1, Fb.kappa(), Fa.kappa());
  // Sum over q = -L, ..., L, where L = 1, so q = -1, 0, 1
  // But, for the Wigner 3j symbol to be non-zero, m1 + m2 + m3 = 0, so the only
  // non-zero term in the sum is for when q = 0.

  double matrix_elements = R*C; //*pow(-1,0.5*(Fa.twoj()-Fb.twoj()))*Angular::threej_2(Fb.twoj(),2,Fa.twoj(),-1,0,1);

  double omega = (Fb.en() - Fa.en());

  // Want just the reduced matrix element squared, remove (1/3)(alpha*omega)^2 term
  return (abs(pow(matrix_elements, 2))*omega / 3.0);
}

double RME_total_e(const HF::HartreeFock *vHF, const DiracSpinor &Fa, double E) {

  const auto &grid = Fa.grid();
  using namespace qip::overloads;
  
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

  // Calculating the form factor of a single state, and then
  // summing over the final states states
  double RME_tot = 0.0;
  for (const auto &Fb : cntm.orbitals) {
    RME_tot += RME_single_e(Fa, Fb, grid.r());
  }
  return RME_tot;
}

// VELOCITY FORM

double RME_single_e_v(const DiracSpinor &Fa, const DiracSpinor &Fb, const auto r) {

  const auto &grid = Fa.grid();
  using namespace qip::overloads;

  std::vector<double> ones(Fa.grid().size(),1);

  // Radial component int (fa*fb + ga*gb) dr
  double V = radial_int(Fb.f(), Fa.g(), ones, grid);
  double W = radial_int(Fb.g(), Fa.f(), ones, grid);

  // Angular component <b|| sigma || a> = 2*<b||s||a>
  double C = 2*Angular::S_kk(-Fb.kappa(), Fa.kappa());
  // Sum over q = -L, ..., L, where L = 1, so q = -1, 0, 1
  // But, for the Wigner 3j symbol to be non-zero, m1 + m2 + m3 = 0, so the only
  // non-zero term in the sum is for when q = 0.

  double matrix_elements = (0.5)*(V*(Fb.kappa()-Fa.kappa()-1) - W*(Fb.kappa()-Fa.kappa()+1)*C);

  double omega = (Fb.en() - Fa.en());

  // Want just the reduced matrix element squared, remove (1/3)(alpha*omega)^2 term
  return (abs(pow(matrix_elements, 2))/ (3.0*omega*omega)); // THIS SHOULD HAVE OMEGA^2 of DENOMINATOR
}

double RME_total_e_v(const HF::HartreeFock *vHF, const DiracSpinor &Fa, double E) {

  const auto &grid = Fa.grid();
  using namespace qip::overloads;
  
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

  // Calculating the form factor of a single state, and then
  // summing over the final states states
  double RME_tot = 0.0;
  for (const auto &Fb : cntm.orbitals) {
    RME_tot += RME_single_e_v(Fa, Fb, grid.r());
  }
  return RME_tot;
}

/*

// OLD ANGULAR INTEGRALS


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

*/
