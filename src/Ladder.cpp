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
void CreationLO(double& ci, int& ni)
{
  //Creation ladder operator
  ci *= sqrt(ni+1); //Update coefficient
  ni += 1; //Update state
  return;
};

void AnnihilationLO(double& ci, int& ni)
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
  for (unsigned int i=0;i<Basis.size();i++)
  {
    //Loop over all modes
    for (int j=0;j<Basis[i].M;j++)
    {
      //Calculate partial energies
      double Ej = 0.5;
      Ej += Basis[i].Modes[j].Quanta;
      Ej *= Basis[i].Modes[j].Freq;
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

