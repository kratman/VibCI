/*

###############################################################################
#                                                                             #
#            Ladder Operator Vibrational Configuration Interaction            #
#                              By: Eric G. Kratz                              #
#                                                                             #
###############################################################################

 Ladder operator routines for LOVCI

*/

//Ladder operators
inline void CreationLO(double& ci, int& ni)
{
  //Creation ladder operator
  ci *= sqrt(ni+1); //Update coefficient
  ni += 1; //Update state
  return;
};

inline void AnnihilationLO(double& ci, int& ni)
{
  //Annihilation ladder operator
  ci *= sqrt(ni); //Update coefficient
  ni -= 1; //Update state
  //Check for impossible states
  if (ni < 0)
  {
    ni = -1; //Used later to remove the state
    ci = 0; //Delete state
  }
  return;
};

//Hamiltonian operators
void ZerothHam(MatrixXd& H)
{
  //Calculate the harmonic Hamiltonian matrix elements
  double Espec = 0; //Spectator mode ZPE energy
  //Spectator modes
  #pragma omp parallel for reduction(+:Espec)
  for (unsigned int i=0;i<SpectModes.size();i++)
  {
    //ZPE for mode i
    double Ei = 0.5*SpectModes[i].Freq;
    //Update energy
    Espec += Ei;
  }
  #pragma omp barrier
  //Harmonic matrix elements
  #pragma omp parallel for
  for (unsigned int i=0;i<BasisSet.size();i++)
  {
    //Loop over all modes
    double Ei = Espec; //Hii matrix element
    for (int j=0;j<BasisSet[i].M;j++)
    {
      //Calculate partial energies
      double Ej = 0.5;
      Ej += BasisSet[i].Modes[j].Quanta;
      Ej *= BasisSet[i].Modes[j].Freq;
      //Update matrix element
      Ei += Ej;
    }
    //Update Hamiltonian
    H(i,i) += Ei;
  }
  #pragma omp barrier
  return;
};

void AnharmHam(MatrixXd& H)
{
  //Add anharmonic terms to the Hamiltonian
  
  return;
};

inline void VCIDiagonalize(MatrixXd& H, MatrixXd& Psi, VectorXd& E)
{
  //Wrapper for the Eigen diagonalization
  EigenSolver<MatrixXd> SE; //Schrodinger equation
  SE.compute(H);
  E = SE.eigenvalues().real();
  Psi = SE.eigenvectors().real();
  return;
};

inline double LBroaden(double fi, double f, double wid)
{
  //Function to calculate the unnormalized Lorentz width
  double lint; //Lorentz intensity
  lint = (fi-f);
  lint *= lint;
  lint += (wid*wid);
  lint = wid/(pi*lint);
  //Return final intensity
  return lint;
};

int IsFund(WaveFunction& bfunc)
{
  //Tests if a basis functionis a fundamental transition
  int fund = -1;
  int sum = 0;
  #pragma omp parallel for reduction(+:sum)
  for (int i=0;i<bfunc.M;i++)
  {
    sum += bfunc.Modes[i].Quanta;
  }
  #pragma omp barrier
  if (sum == 1)
  {
    #pragma omp parallel for
    for (int i=0;i<bfunc.M;i++)
    {
      if (bfunc.Modes[i].Quanta == 1)
      {
        fund = i;
      }
    }
    #pragma omp barrier
  }
  return fund;
};

void PrintSpectrum(VectorXd& Freqs, fstream& outfile)
{
  //Function to print the CI spectrum
  double Fmin = 0; //Start of the spectrum
  double Fmax = Freqs.maxCoeff()+(20*LorentzWid); //End of the spectrum
  double fn = Fmin; //Current frequency at point n
  //All intensity derives from the fundamental transitions
  while (fn <= Fmax)
  {
    double In = 0; //Intensity at point n
    //Spectator modes
    #pragma omp parallel for reduction(+:In)
    for (unsigned int i=0;i<SpectModes.size();i++)
    {
      double I;
      //All spectators are fundamentals
      I = SpectModes[i].ModeInt;
      I *= LBroaden(fn,SpectModes[i].Freq,LorentzWid);
      //Sum intensities
      In += I;
    }
    #pragma omp barrier
    //VCI modes
    #pragma omp parallel for reduction(+:In)
    for (unsigned int i=0;i<BasisSet.size();i++)
    {
      //Add fundamental intensities
      if (Freqs(i) > 0)
      {
        //Add all modes besides the VCI ground state
        
      }
    }
    #pragma omp barrier
    //Write data
    outfile << fn << " " << In << '\n';
    //Go to point n+1
    fn += DeltaFreq;
  }
  return;
};

