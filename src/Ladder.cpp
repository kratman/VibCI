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

double AnharmPot(int n, int m, FConst& fc)
{
  //Calculate anharmonic matrix elements for <m|H|n>
  double Vnm = 0;
  vector<vector<int> > NewStates;
  vector<vector<double> > StateCoeffs;
  //Sort force constants
  vector<int> ShortModes; //A short list of the key modes
  vector<int> ModePowers; //Powers in each mode
  for (unsigned int i=0;i<fc.fcpow.size();i++)
  {
    //Sort the modes and powers
    if (i == 0)
    {
      ShortModes.push_back(fc.fcpow[i]);
      ModePowers.push_back(1);
      vector<int> tmp1;
      vector<double> tmp2;
      NewStates.push_back(tmp1);
      StateCoeffs.push_back(tmp2);
    }
    else
    {
      bool foundmode = 0;
      for (unsigned int j=0;j<ShortModes.size();j++)
      {
        if (fc.fcpow[i] == ShortModes[j])
        {
          ModePowers[j] += 1;
          foundmode = 1;
        }
      }
      if (!foundmode)
      {
        ShortModes.push_back(fc.fcpow[i]);
        ModePowers.push_back(1);
        vector<int> tmp1;
        vector<double> tmp2;
        NewStates.push_back(tmp1);
        StateCoeffs.push_back(tmp2);
      }
    }
  }
  //Create new states
  for (unsigned int i=0;i<ShortModes.size();i++)
  {
    //Apply operators
    NewStates[i].push_back(BasisSet[n].Modes[ShortModes[i]].Quanta);
    StateCoeffs[i].push_back(1.0);
    for (int j=0;j<ModePowers[i];j++)
    {
      //Create new state for mode i
      vector<int> stateupdate;
      vector<double> coeffupdate;
      for (unsigned int k=0;k<NewStates[i].size();k++)
      {
        int quant;
        double coeff;
        //Creation
        quant = NewStates[i][k];
        coeff = StateCoeffs[i][k];
        CreationLO(coeff,quant);
        stateupdate.push_back(quant);
        coeffupdate.push_back(coeff);
        //Annihilation
        quant = NewStates[i][k];
        coeff = StateCoeffs[i][k];
        AnnihilationLO(coeff,quant);
        if (quant >= 0)
        {
          stateupdate.push_back(quant);
          coeffupdate.push_back(coeff);
        }
      }
      //Save states
      NewStates[i] = stateupdate;
      StateCoeffs[i] = coeffupdate;
    }
  }
  //Sum energies
  vector<double> CoeffSum;
  for (unsigned int i=0;i<ShortModes.size();i++)
  {
    CoeffSum.push_back(0.0);
    for (unsigned int j=0;j<NewStates[i].size();j++)
    {
      int quantn = NewStates[i][j];
      int quantm = BasisSet[m].Modes[ShortModes[i]].Quanta;
      if (quantn == quantm)
      {
        CoeffSum[i] += StateCoeffs[i][j];
      }
    }
  }
  //Scale by the force constant
  Vnm = fc.fc;
  //Combine coeffcients
  for (unsigned int i=0;i<CoeffSum.size();i++)
  {
    Vnm *= CoeffSum[i];
  }
  return Vnm;
};

//Hamiltonian operators
void ZerothHam(MatrixXd& H)
{
  //Calculate the harmonic Hamiltonian matrix elements
  double Espec = 0; //Spectator mode ZPE energy
  //Spectator modes
  for (unsigned int i=0;i<SpectModes.size();i++)
  {
    //ZPE for mode i
    double Ei = 0.5*SpectModes[i].Freq;
    //Update energy
    Espec += Ei;
  }
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
  #pragma omp parallel for
  for (unsigned int i=0;i<BasisSet.size();i++)
  {
    for (unsigned int j=0;j<BasisSet.size();j++)
    {
      double Vij = 0;
      for (unsigned int k=0;k<AnharmFC.size();k++)
      {
        if (ScreenState(i,j,AnharmFC[k]))
        {
          //Add anharmonic matrix elements
          Vij += AnharmPot(i,j,AnharmFC[k]);
        }
      }
      H(i,j) += Vij;
    }
  }
  #pragma omp barrier
  return;
};

//Utility functions
inline void VCIDiagonalize(MatrixXd& H, MatrixXd& Psi, VectorXd& E)
{
  //Wrapper for the Eigen diagonalization
  EigenSolver<MatrixXd> SE; //Schrodinger equation
  SE.compute(H); //Diagonalize the matrix
  E = SE.eigenvalues().real(); //Extract frequencies
  Psi = SE.eigenvectors().real(); //Extract CI vectors
  return;
};

inline double LBroaden(double fi, double f, double wid)
{
  //Function to calculate the Lorentz width
  double lint; //Lorentz intensity
  lint = (fi-f);
  lint *= lint;
  lint += (wid*wid);
  lint = wid/(pi*lint);
  //Return final intensity
  return lint;
};

inline double GBroaden(double fi, double f, double wid)
{
  //Function to calculate the Gaussian width
  double gint; //Lorentz intensity
  gint = (fi-f);
  gint *= gint;
  gint /= (2*wid*wid);
  gint = exp(-1*gint);
  gint /= (wid*rt2pi);
  //Return final intensity
  return gint;
};

int IsFund(WaveFunction& bfunc)
{
  //Tests if a basis functionis a fundamental transition
  int fund = -1;
  int sum = 0;
  for (int i=0;i<bfunc.M;i++)
  {
    sum += bfunc.Modes[i].Quanta;
  }
  if (sum == 1)
  {
    for (int i=0;i<bfunc.M;i++)
    {
      if (bfunc.Modes[i].Quanta == 1)
      {
        fund = i;
      }
    }
  }
  return fund;
};

bool ScreenState(int n, int m, FConst& fc)
{
  //Function for ignoring states that have no overlap
  bool keepstate = 1; //Assume the state is good
  //Sort force constants
  vector<int> ShortModes; //A short list of the key modes
  vector<int> ModePowers; //Powers in each mode
  for (unsigned int i=0;i<fc.fcpow.size();i++)
  {
    //Sort the modes and powers
    if (i == 0)
    {
      //Initialize the vector
      ShortModes.push_back(fc.fcpow[i]);
      ModePowers.push_back(1);
    }
    else
    {
      //Update powers
      bool foundmode = 0;
      for (unsigned int j=0;j<ShortModes.size();j++)
      {
        if (fc.fcpow[i] == ShortModes[j])
        {
          ModePowers[j] += 1;
          foundmode = 1;
        }
      }
      if (!foundmode)
      {
        ShortModes.push_back(fc.fcpow[i]);
        ModePowers.push_back(1);
      }
    }
  }
  //Check based on force constant powers
  for (unsigned int i=0;i<ShortModes.size();i++)
  {
    int qdiff = 0; //Number of quanta between states
    qdiff += BasisSet[n].Modes[ShortModes[i]].Quanta;
    qdiff -= BasisSet[m].Modes[ShortModes[i]].Quanta;
    qdiff = abs(qdiff);
    if (qdiff > ModePowers[i])
    {
      //Impossible for the states to overlap
      keepstate = 0;
    }
  }
  if (!keepstate)
  {
    //Skip the rest of the checks
    return 0;
  }
  //Check overlap of all other modes
  for (int i=0;i<BasisSet[n].M;i++)
  {
    bool cont = 1; //Continue the check
    for (unsigned int j=0;j<ShortModes.size();j++)
    {
      if (ShortModes[j] == i)
      {
        //Ignore this mode since it is in the FC
        cont = 0;
      }
    }
    if (cont)
    {
      if (BasisSet[n].Modes[i].Quanta != BasisSet[m].Modes[i].Quanta)
      {
        //Remove state due to zero overlap in mode i
        keepstate = 0;
      }
    }
  }
  //Return decision
  return keepstate;
};

void PrintSpectrum(VectorXd& Freqs, MatrixXd& Psi, fstream& outfile)
{
  //Function to print the CI spectrum
  double Fmin = 0; //Start of the spectrum
  double Fmax = Freqs.maxCoeff()+(20*LorentzWid); //End of the spectrum
  if (Fmax > FreqCut)
  {
    //Scale maximum frequency
    Fmax = FreqCut;
  }
  double fn = Fmin; //Current frequency at point n
  //Look for threading errors
  int SpectThreads = Ncpus;
  if (SpectModes.size() <= 1)
  {
    //Run serial
    SpectThreads = 1;
  }
  //All intensity derives from the fundamental transitions
  while (fn <= Fmax)
  {
    double In = 0; //Intensity at point n
    //Spectator modes
    #pragma omp parallel for reduction(+:In) num_threads(SpectThreads)
    for (unsigned int i=0;i<SpectModes.size();i++)
    {
      double I;
      //All spectators are fundamentals
      I = SpectModes[i].ModeInt;
      if (GauBroad)
      {
        I *= GBroaden(fn,SpectModes[i].Freq,LorentzWid);
      }
      else
      {
        I *= LBroaden(fn,SpectModes[i].Freq,LorentzWid);
      }
      //Sum intensities
      In += I;
    }
    #pragma omp barrier
    //VCI modes
    for (unsigned int i=0;i<BasisSet.size();i++)
    {
      //Loop over states
      if ((Freqs(i) > 0) and (Freqs(i) < (1.10*Fmax)))
      {
        //Add all states besides the VCI ground state
        #pragma omp parallel for reduction(+:In)
        for (unsigned int j=0;j<BasisSet.size();j++)
        {
          //Loop over basis functions and find fundamentals
          int FundID = IsFund(BasisSet[j]);
          if (FundID >= 0)
          {
            double I;
            //Scale by intensity
            I = BasisSet[j].Modes[FundID].ModeInt;
            //Add weight
            I *= Psi(j,i)*Psi(j,i);
            //Broaden frequency
            if (GauBroad)
            {
              I *= GBroaden(fn,Freqs(i),LorentzWid);
            }
            else
            {
              I *= LBroaden(fn,Freqs(i),LorentzWid);
            }
            //Sum intensities
            In += I;
          }
        }
        #pragma omp barrier
      }
    }
    //Write data
    outfile << fn << " " << In << '\n';
    //Go to point n+1
    fn += DeltaFreq;
  }
  //Print and return
  outfile.flush();
  return;
};

void ScaleFC()
{
  //Adjust force constants for permutations and sqrt(2) terms
  for (unsigned int i=0;i<AnharmFC.size();i++)
  {
    //Modify all force constants
    double fcscale;
    //Add sqrt(2) terms
    fcscale = ((double)AnharmFC[i].fcpow.size());
    fcscale = pow(2,fcscale);
    fcscale = 1/sqrt(fcscale);
    //Sort force constants
    vector<int> ShortModes; //A short list of the key modes
    vector<int> ModePowers; //Powers in each mode
    for (unsigned int j=0;j<AnharmFC[i].fcpow.size();j++)
    {
      //Sort the modes and powers
      if (j == 0)
      {
        //Initialize the vector
        ShortModes.push_back(AnharmFC[i].fcpow[j]);
        ModePowers.push_back(1);
      }
      else
      {
        //Update powers
        bool foundmode = 0;
        for (unsigned int k=0;k<ShortModes.size();k++)
        {
          if (AnharmFC[i].fcpow[j] == ShortModes[k])
          {
            ModePowers[k] += 1;
            foundmode = 1;
          }
        }
        if (!foundmode)
        {
          ShortModes.push_back(AnharmFC[i].fcpow[j]);
          ModePowers.push_back(1);
        }
      }
    }
    //Calculate permutations
    for (unsigned int j=0;j<ModePowers.size();j++)
    {
      fcscale /= Fact(ModePowers[j]);
    }
    //Scale FC
    AnharmFC[i].fc *= fcscale;
  }
  return;
};

