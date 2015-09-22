/*

###############################################################################
#                                                                             #
#            Ladder Operator Vibrational Configuration Interaction            #
#                              By: Eric G. Kratz                              #
#                                                                             #
###############################################################################

 Headers, libraries, and data structures for LOVCI

*/

//Make including safe
#ifndef VCI_HEADERS
#define VCI_HEADERS

//Header Files
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <complex>
#include <cmath>
#include <fstream>
#include <vector>
#include <map>
#include <sys/stat.h>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <Eigen/StdList>
#include <Eigen/Eigen>
#include <Eigen/StdVector>

//Set namespaces for common libraries
using namespace Eigen;
using namespace std;

//Global exact constants
const double pi = 4*atan(1); //Pi

//Global measured constants (NIST, CODATA 2010)
const double cs = 2.99792458e-10; //Speed of light (cm)
const double k = 0.69503476; //Boltzmann constant (cm^-1)

//Global derived constants
const double h = 2*pi; //Planck constant (cm^-1)

//Timers
int StartTime = 0; //Time the calculation starts
int EndTime = 0; //Time the calculation ends

//Custom classes
class HOFunc
{
  //Class for harmonic oscillator basis functions
  public:
    bool Spect; //Spectator mode
    double Freq; //Frequency
    int Quanta; //Number of quanta in the mode
    double ModeInt; //Intensity
};

class WaveFunction
{
  //Class for storing a VCI wavefunction
  public:
    int M; //Number of modes
    vector<HOFunc> Modes; //Functions
};

//Global variables
int Ncpus = 0; //Number of CPUs for the calculations
double LorentzWid = 1; //Width of the peaks in the final spectrum
vector<WaveFunction> BasisSet; //Full basis set

//Function declarations (alphabetical)
void AnharmHam(MatrixXd&);

void AnnihilationLO(double&,int&);

bool CheckFile(const string&);

void CreationLO(double&,int&);

int FindMaxThreads();

double LBroaden(double,double,double);

void PrintFancyTitle();

void ReadCIArgs(int,char*,fstream&,fstream&);

void ReadCIInput(MatrixXd&,fstream&);

void ZerothHam(MatrixXd&);

//Function definitions (alphabetical)
#include "Core_functions.cpp"
#include "Input_Reader.cpp"
#include "Ladder.cpp"

#endif
