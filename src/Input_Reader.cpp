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
  bool DoQuit = 0; //Exit if an error is detected
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
  //Check for argument errors
  if (Ncpus < 1)
  {
    //Checks the number of threads and continue
    cout << " Warning: Calculations cannot run with ";
    cout << Ncpus << " CPUs.";
    cout << '\n';
    cout << "  Do you know how computers work?";
    cout << " Ncpus set to 1";
    cout << '\n';
    Ncpus = 1;
    cout.flush(); //Print warning
  }
  if (Ncpus > FindMaxThreads())
  {
    cout << " Warning: Too many threads requested.";
    cout << '\n';
    cout << "  Requested: " << Ncpus << ",";
    cout << " Available: " << FindMaxThreads();
    cout << '\n';
    cout << "  Ncpus set to " << FindMaxThreads();
    cout << '\n';
    Ncpus = FindMaxThreads();
    cout.flush(); //Print warning
  }
  if (!vcidata.good())
  {
    //Check input file
    cout << " Error: Could not open input file.";
    cout << '\n';
    DoQuit = 1;
  }
  if (!spectfile.good())
  {
    //Check output file
    cout << " Error: Could not output file.";
    cout << '\n';
    DoQuit = 1;
  }

  if (DoQuit)
  {
    //Quit if there is an error
    cout << '\n';
    cout.flush();
    exit(0);
  }
  else
  {
    //Sarcastically continue
    cout << '\n';
    cout << "No fatal errors detected.";
    cout << '\n';
    cout << " And there was much rejoicing. Yay...";
    cout << '\n' << '\n';
    cout.flush();
  }
  //Set threads
  omp_set_num_threads(Ncpus);
  return;
};

void ReadCIInput(MatrixXd& VCIHam, fstream& vcidata)
{
  //Function to read the input files
  string dummy; //Generic sting
  //Count basis functions and read modes
  bool ProgSet = 0; //Make the basis a progression in a mode
  unsigned int progmode = 0; //Mode for the progression
  int ProgQuanta = 0; //Number of quanta for the progression
  vector<unsigned int> ProgModes; //Modes with progressions
  int Nspect = 0; //Number of spectator modes
  int Nmodes = 0; //Number of different modes
  vector<HOFunc> BasisCount; //Temp. storage of modes
  int Nfc = 0; //Number of force constants
  //Read basis set type
  vcidata >> dummy >> dummy;
  if ((dummy == "Progression") or (dummy == "progression"))
  {
    ProgSet = 1;
    vcidata >> dummy; //Clear junk
    vcidata >> progmode; //Mode that forms the progressions
    vcidata >> ProgQuanta; //Number of quanta in the progression
    int modect; //Number of modes with progressions
    vcidata >> dummy >> modect;
    for (int i=0;i<modect;i++)
    {
      int tmp;
      vcidata >> tmp;
      ProgModes.push_back(tmp);
    }
  }
  //Read broadening settings
  vcidata >> dummy >> LorentzWid >> DeltaFreq;
  //Read active modes
  vcidata >> dummy; //Clear junk
  vcidata >> Nmodes; //Read modes
  for (int i=0;i<Nmodes;i++)
  {
    //Actual vibrational modes
    HOFunc tmp;
    vcidata >> dummy; //Throw out ID number
    vcidata >> tmp.Freq; //Frequency
    vcidata >> tmp.Quanta; //Max number of quanta
    vcidata >> tmp.ModeInt; //Intensity
    BasisCount.push_back(tmp);
  }
  //Read spectator modes
  vcidata >> dummy; //Clear junk
  vcidata >> Nspect; //Read spectator modes
  for (int i=0;i<Nspect;i++)
  {
    //Spectator modes
    HOFunc tmp;
    vcidata >> dummy; //Throw out ID number
    vcidata >> tmp.Freq; //Frequency
    tmp.Quanta = 1; //Only one state
    vcidata >> tmp.ModeInt; //Intensity
    SpectModes.push_back(tmp);
  }
  //Count states
  if (ProgSet)
  {
    //Simple progression basis
    Nmodes = 1;
    for (unsigned int i=0;i<BasisCount.size();i++)
    {
      Nmodes += BasisCount[i].Quanta;
    }
    Nmodes += ProgQuanta*ProgModes.size();
  }
  else
  {
    //Product basis
    for (unsigned int i=0;i<BasisCount.size();i++)
    {
      if (i == 0)
      {
        Nmodes = (BasisCount[i].Quanta+1);
      }
      else
      {
        Nmodes *= (BasisCount[i].Quanta+1);
      }
    }
  }
  //Read anharmonic force constants
  vcidata >> dummy; //Clear junk
  vcidata >> Nfc; //Read number of anharmonic force constants
  for (int i=0;i<Nfc;i++)
  {
    //Save force constant data
    FConst tmp;
    int fcpower = 0;
    vcidata >> fcpower;
    for (int j=0;j<fcpower;j++)
    {
      int modej = 0;
      vcidata >> modej;
      tmp.fcpow.push_back(modej);
    }
    vcidata >> tmp.fc; //Read force constant value
    AnharmFC.push_back(tmp);
  }
  //Create data structures
  VCIHam = MatrixXd(Nmodes,Nmodes); //Create the Hamiltonian matrix
  VCIHam.setZero(); //Ensure that the Hamiltonian matrix is empty
  if (ProgSet)
  {
    //Create progression basis
    WaveFunction GroundState;
    GroundState.M = BasisCount.size();
    for (unsigned int k=0;k<BasisCount.size();k++)
    {
      //Copy general information
      HOFunc tmp;
      tmp.Freq = BasisCount[k].Freq;
      tmp.ModeInt = BasisCount[k].ModeInt;
      tmp.Quanta = 0;
      //Save basis function component
      GroundState.Modes.push_back(tmp);
    }
    BasisSet.push_back(GroundState); //Ground state reference
    for (unsigned int i=0;i<BasisCount.size();i++)
    {
      //Add 1D modes
      for (int j=0;j<BasisCount[i].Quanta;j++)
      {
        WaveFunction temp;
        temp.M = BasisCount.size();
        for (unsigned int k=0;k<BasisCount.size();k++)
        {
          //Copy general information
          HOFunc tmp;
          tmp.Freq = BasisCount[k].Freq;
          tmp.ModeInt = BasisCount[k].ModeInt;
          //Add quanta
          if (i == k)
          {
            tmp.Quanta = (j+1);
          }
          else
          {
            tmp.Quanta = 0;
          }
          //Save basis function component
          temp.Modes.push_back(tmp);
        }
      BasisSet.push_back(temp);
      }
    }
    //Add combination bands
    for (unsigned int i=0;i<ProgModes.size();i++)
    {
      //Add 1D modes
      for (int j=0;j<ProgQuanta;j++)
      {
        WaveFunction temp;
        temp.M = BasisCount.size();
        for (unsigned int k=0;k<BasisCount.size();k++)
        {
          //Copy general information
          HOFunc tmp;
          tmp.Freq = BasisCount[k].Freq;
          tmp.ModeInt = BasisCount[k].ModeInt;
          //Add quanta
          if (i == progmode)
          {
            tmp.Quanta = (j+1);
          }
          else if (i == ProgModes[i])
          {
            tmp.Quanta = 1;
          }
          else
          {
            tmp.Quanta = 0;
          }
          //Save basis function component
          temp.Modes.push_back(tmp);
        }
      BasisSet.push_back(temp);
      }
    }
  }
  else
  {
    //Create product basis
    for (int i=0;i<(BasisCount[0].Quanta+1);i++)
    {
      WaveFunction temp;
      HOFunc tmp;
      tmp.Freq = BasisCount[0].Freq;
      tmp.ModeInt = BasisCount[0].ModeInt;
      tmp.Quanta = i;
      temp.Modes.push_back(tmp);
      BasisSet.push_back(temp);
    }
    for (unsigned int i=1;i<BasisCount.size();i++)
    {
      vector<WaveFunction> NewBasis;
      for (int j=0;j<(BasisCount[i].Quanta+1);j++)
      {
        vector<WaveFunction> BasisCopy = BasisSet;
        for (unsigned int k=0;k<BasisSet.size();k++)
        {
          HOFunc tmp;
          //Copy data
          tmp.Freq = BasisCount[i].Freq;
          tmp.ModeInt = BasisCount[i].ModeInt;
          //Set quanta
          tmp.Quanta = j;
          //Add to arrays
          BasisCopy[k].Modes.push_back(tmp);
        }
        for (unsigned int n=0;n<BasisCopy.size();n++)
        {
          NewBasis.push_back(BasisCopy[n]);
        }
      }
      BasisSet = NewBasis;
    }
  }
  //Correct array lengths
  #pragma omp parallel for
  for (unsigned int i=0;i<BasisSet.size();i++)
  {
    //Update counter
    BasisSet[i].M = BasisSet[i].Modes.size();
  }
  #pragma omp barrier
  //Print settings
  cout << "General settings:" << '\n';
  cout << "  Threads: " << Ncpus << '\n';
  cout << '\n';
  cout << "Spectrum settings:" << '\n';
  cout << "  Active modes: " << BasisCount.size() << '\n';
  cout << "  Spectator modes: " << SpectModes.size() << '\n';
  cout.precision(3); //Truncate numbers
  cout << "  Line width: " << LorentzWid << '\n';
  cout << "  Resolution: " << DeltaFreq << '\n';
  cout.precision(12); //Replace settings
  cout << '\n';
  cout << "VCI basis set:" << '\n';
  cout << "  Basis type: ";
  if (ProgSet)
  {
    cout << "Progression" << '\n';
    cout << "    Mode ID: " << progmode << '\n';
    cout.precision(2); //Truncate frequency
    cout << "    Frequency: " << BasisCount[progmode].Freq << '\n';
    cout.precision(12); //Replace settings
    cout << "    Quanta: " << ProgQuanta << '\n';
  }
  else
  {
    cout << "Product" << '\n';
  }
  cout << "  Basis functions: " << BasisSet.size() << '\n';
  cout << "  Force constants: ";
  if (Nfc == 0)
  {
    cout << "N/A" << '\n';
  }
  else
  {
    cout << AnharmFC.size() << '\n';
  }
  cout << '\n';
  cout.flush();
  return;
};

