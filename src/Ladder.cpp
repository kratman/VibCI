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
  #pragma omp parallel for
  for (unsigned int i=0;i<BasisSet.size();i++)
  {
    //Loop over all modes
    for (int j=0;j<BasisSet[i].M;j++)
    {
      //Calculate partial energies
      double Ej = 0.5;
      Ej += BasisSet[i].Modes[j].Quanta;
      Ej *= BasisSet[i].Modes[j].Freq;
      //Update Hamiltonian
      H(i,i) += Ej;
    }
  }
  #pragma omp barrier
  return;
};

void AnharmHam(MatrixXd& H)
{
  //Add anharmonic terms to the Hamiltonian
  
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
  return lint;
};
