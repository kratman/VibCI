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
const double sqrt2 = pow(2,0.5); //Square root of 2

//Global measured constants (NIST, CODATA 2010)
const double k = 0; //Boltzmann constant (cm^-1)
const double Har2eV = 27.21138505; //Hartrees to eV

//Globals


//Timers
int StartTime = 0; //Time the calculation starts
int EndTime = 0; //Time the calculation ends

//Custom classes


//Function declarations (alphabetical)


//Function definitions (alphabetical)


#endif