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
  //End of section

  //Initialize local variables
  fstream vcidata,spectfile; //File streams for the input and output files
  MatrixXd VCIHam; //Full Hamiltonian matrix
  double RunTime; //Run time for the calculations
  string TimeUnits; //Seconds, minutes, or hours for RunTime
  //End of section

  //Print title and compile date
  PrintFancyTitle();
  cout << '\n';
  cout << "Last modification: ";
  cout << __TIME__ << " on ";
  cout << __DATE__ << '\n';
  cout << '\n';
  cout.flush();
  //End of section

  //Gather input and check for errors
  ReadCIArgs(argc,argv,vcidata,spectfile); //Read arguments
  ReadCIInput(VCIHam,vcidata); //Read input files
  //End of section

  //Calculate spectrum
  ZerothHam(VCIHam);
  AnharmHam(VCIHam);
  //End of section

  //Print results
  cout << '\n';
  
  EndTime = (unsigned)time(0); //Time the program starts
  RunTime = (double)(EndTime-StartTime);
  if (RunTime >= 3600)
  {
    RunTime /= 3600;
    TimeUnits = "hours";
  }
  else if (RunTime >= 60)
  {
    RunTime /= 60;
    TimeUnits = "minutes";
  }
  else
  {
    TimeUnits = "seconds";
  }
  cout << "Run time: " << RunTime;
  cout << " " << TimeUnits;
  cout << '\n' << '\n';
  cout << "Done.";
  cout << '\n' << '\n';
  //End of section

  //Quit
  return 0;
};
