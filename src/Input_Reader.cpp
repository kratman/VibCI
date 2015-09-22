/*

###############################################################################
#                                                                             #
#            Ladder Operator Vibrational Configuration Interaction            #
#                              By: Eric G. Kratz                              #
#                                                                             #
###############################################################################

 Input reader and error checker for LOVCI

*/

//Input reader
void ReadCIArgs(int argc, char* argv[], fstream& vcidata, fstream& spectfile)
{
  //Function to read the command line arguments
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  //Read command line arguments
  if (argc == 1)
  {
    //Escape if there are no arguments
    cout << '\n';
    cout << "Missing arguments...";
    cout << '\n' << '\n';
    cout << "Usage: lovci -n Ncpus -i Modes.inp -o Spectrum.txt ";
    cout << '\n' << '\n';
    cout << "Use -h or --help for detailed instructions.";
    cout << '\n' << '\n';
    cout.flush();
    exit(0);
  }
  if ((argc % 2) != 1)
  {
    dummy = string(argv[1]);
    if ((dummy != "-h") and (dummy != "--help"))
    {
      //Escape if there are missing arguments
      cout << '\n';
      cout << "Wrong number of arguments...";
      cout << '\n' << '\n';
      cout << "Usage: lovci -n Ncpus -i Modes.inp -o Spectrum.txt ";
      cout << '\n' << '\n';
      cout << "Use -h or --help for detailed instructions.";
      cout << '\n' << '\n';
      cout.flush();
      exit(0);
    }
  }
  for (int i=0;i<argc;i++)
  {
    //Read file names and CPUs
    dummy = string(argv[i]);
    if ((dummy == "-h") or (dummy == "--help"))
    {
      //Print helpful information and exit
      cout << '\n';
      cout << "Usage: lovci -n Ncpus -i Modes.inp -o Spectrum.txt ";
      cout << '\n' << '\n';
      cout << "Command line arguments:";
      cout << '\n' << '\n';
      cout << "  -n    Number of CPUs to use for the calculation.";
      cout << '\n' << '\n';
      cout << "  -i    Input file.";
      cout << '\n' << '\n';
      cout << "  -o    Output file for the VCI spectrum.";
      cout << '\n' << '\n';
      cout.flush();
      exit(0);
    }
    if (dummy == "-n")
    {
      Ncpus = atoi(argv[i+1]);
    }
    if (dummy == "-i")
    {
      vcidata.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-o")
    {
      spectfile.open(argv[i+1],ios_base::out);
    }
  }
  return;
};

void ReadCIInput(MatrixXd& VCIHam, fstream& vcidata)
{
  //Function to read the input files
  
  return;
};

