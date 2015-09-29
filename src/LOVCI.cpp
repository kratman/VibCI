/*

###############################################################################
#                                                                             #
#            Ladder Operator Vibrational Configuration Interaction            #
#                              By: Eric G. Kratz                              #
#                                                                             #
###############################################################################

 LOVCI is licensed under GPLv3, for more information see GPL_LICENSE

*/

//Main header
#include "VCI_headers.h"

int main(int argc, char* argv[])
{
  //Misc. initialization
  StartTime = (unsigned)time(0); //Time the program starts
  cout.precision(12);
  cout << fixed;
  //End of section

  //Initialize local variables
  fstream vcidata,spectfile; //File streams for the input and output files
  MatrixXd VCIHam; //Full Hamiltonian matrix
  MatrixXd CIVec; //CI eigenvectors (wavefunction)
  VectorXd CIFreq; //CI eigenvalues (frequencies
  double Ezpe = 0; //CI zero-point energy
  double RunTime; //Run time for the calculations
  string TimeUnits; //Seconds, minutes, or hours for RunTime
  //End of section

  //Print title and compile date
  PrintFancyTitle();
  cout << "Last modification: ";
  cout << __TIME__ << " on ";
  cout << __DATE__ << '\n';
  cout << '\n';
  cout.flush();
  //End of section

  //Gather input and check for errors
  cout << "Reading input..." << '\n';
  ReadCIArgs(argc,argv,vcidata,spectfile); //Read arguments
  ReadCIInput(VCIHam,vcidata); //Read input files
  //End of section

  //Calculate spectrum
  cout << "Constructing the CI Hamiltonian...";
  cout << '\n';
  cout.flush(); //Print progress
  cout << "  Zeroth-order Hamiltionian";
  ZerothHam(VCIHam); //Harmonic terms
  cout << "; Done." << '\n';
  cout << "  Adding anharmonic potential";
  AnharmHam(VCIHam); //Anharmonic terms
  cout << "; Done." << '\n';
  cout << '\n';
  cout << "Diagonalizing the Hamiltonian...";
  cout << '\n' << '\n';
  cout.flush(); //Print progress
  VCIDiagonalize(VCIHam,CIVec,CIFreq); //Diagonalization wrapper
  //End of section

  //Calculate and remove ZPE
  VectorXd UnitVec(BasisSet.size());
  UnitVec.setOnes(); //Create a unit vector
  Ezpe = CIFreq.minCoeff(); //Save ZPE
  CIFreq -= (Ezpe*UnitVec); //Remove ZPE
  //End of section

  //Print results
  cout << "Printing the spectrum...";
  cout << '\n' << '\n';
  cout << "Results:" << '\n';
  cout << "  Zero-point energy: ";
  cout.precision(2); //Truncate energy
  cout << Ezpe << " (1/cm)"; //Print energy
  cout.precision(12); //Replace settings
  cout << '\n';
  cout.flush(); //Print progress
  PrintSpectrum(CIFreq,CIVec,spectfile); //Convert the frequencies to a spectrum
  EndTime = (unsigned)time(0); //Time the calculation stops
  RunTime = (double)(EndTime-StartTime); //Total run time
  if (RunTime >= 3600)
  {
    //Switch to hours
    RunTime /= 3600;
    TimeUnits = "hours";
  }
  else if (RunTime >= 60)
  {
    //Switch to minutes
    RunTime /= 60;
    TimeUnits = "minutes";
  }
  else
  {
    //Stick with seconds
    TimeUnits = "seconds";
  }
  cout.precision(2); //Truncate time
  cout << "  Run time: " << RunTime;
  cout.precision(12); //Replace settings
  cout << " " << TimeUnits;
  cout << '\n' << '\n';
  cout << "Done.";
  cout << '\n' << '\n';
  //End of section

  //Quit
  return 0;
};
